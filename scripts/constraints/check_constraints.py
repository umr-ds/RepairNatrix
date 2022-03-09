from scripts.constraints.gc_content import overall_gc_content_error_val
from scripts.constraints.homopolymers import homopolymer_error_val
from scripts.constraints.kmer import kmer_counting_error_val
from scripts.constraints.undesired_subsequences import undesired_subsequences_val
import numpy as np

allowed_bases = {"A", "C", "G", "T"}
DEBUG = True
quality_score_format_dict = {"Illumina 1.8+": 33, "PacBio": 33, "Sanger": 33,
                             "Solexa": 64, "Illumina 1.3+": 64, "Illumina 1.5+": 64}

"""
required config entries:
    - allowed min / max gc content ( <= / >= )
    - allowed max homopolymer length ( <= )
    - inplace repair: if corrupt sequences should remain in the output or be replaced with the repaired sequence  (bool)
    - repair quality-score: phred quality score to set for a changed/inserted/removed base (int)
        a deletion will change the base infront of the deletion to the quality score...
    kmer_k: kmer length (int)
    kmer_max_count: maximum number of kmer occurrences ( <= )
    
"""

inplace_repair = snakemake.params.inplace_repair
repair_quality_score = snakemake.params.repair_quality_score
allowed_min_gc, allowed_max_gc = snakemake.params.allowed_gc
allowed_max_homopolymer_length = snakemake.params.allowed_max_homopolymer_length
kmer_k = snakemake.params.kmer_k
kmer_max_count = snakemake.params.kmer_max_count

""" 
idea:
    - AFTER clustering: check if there are multiple clusters with the same seed (i.e. the same X starting bases)
        - parameter: unique first X bases (int)
    - BEFORE forward/reverse-read-merging: limit GC content rule to windowed mode and skip last/first window.
"""

def allmax(a):
    if len(a) == 0:
        return []
    all_ = [0]
    max_ = a[0]
    for i in range(1, len(a)):
        if a[i] > max_:
            all_ = [i]
            max_ = a[i]
        elif a[i] == max_:
            all_.append(i)
    return all_


def calc_errors(seq):
    """
    Calculates the error-values for each constraint.
    :param seq: sequence to test
    :return: list of lists of error-values for each constraint:
        [undesired_subsequences, overall_gc_content, homopolymer_error, kmer_counting, list of 0 for minus optimization]
    """
    results = [undesired_subsequences_val(seq, "undesired_subsequences.txt"),
               overall_gc_content_error_val(seq, float(allowed_min_gc), float(allowed_max_gc)),
               homopolymer_error_val(seq, int(allowed_max_homopolymer_length), True),
               kmer_counting_error_val(seq, int(kmer_k), int(kmer_max_count)),
               [0.0] * len(seq)]
    return results


def repairInsertionWithDeletion(seq, desired_length):
    # the sequence is a result of (multiple) insertions and deletions
    # while this should not be the first option to choose (orkhams razor), this CAN happen...

    # first we enforce the desired length just like try_repair() does

    # then we try removing a base to lower the error value
    # -> if there is a different location with an error remaining, we then try to insert a base at that location
    # -> if there is no error remaining:
    #   -> the wrong base is removed
    #   -> the missing insertion did not create a rule-violation
    pass


def try_repair(seq, desired_length):
    org_seq = seq
    iterations = 0
    no_changes = 0
    results = calc_errors(seq)
    res = [sum(x) for x in zip(*results)]

    while len(seq) < desired_length:
        # DELETION ERROR:
        # this is a bit more complicated: we cant simply pick the position with the highest error-value
        # because errors might appear _after_ a possible deletion
        # however we _can_ be sure that the deletion should be _before_

        # to/do: right now we randomly take the argmax, but we might want to go for the middle-position (minimizes the number of changes for homopolymers)
        # max_arg = np.random.choice(allmax(res))
        max_arg = allmax(res)
        max_arg = max_arg[int(len(max_arg) / 2)]

        possible_bases = allowed_bases - {seq[max_arg]}
        min_seq = seq
        min_sum_res = sum(res)
        insert = 0
        for base in possible_bases:
            seq = seq[:max_arg] + base + seq[max_arg + insert:]  # for the first run we insert instead of replacing..
            insert = 1
            results = calc_errors(seq)
            res = [sum(x) for x in zip(*results)]
            if sum(res) < min_sum_res:
                min_seq = seq
                min_sum_res = sum(res)
            if sum(res) == 0:
                break
        if seq != min_seq:
            no_changes += 1
            seq = min_seq
            if DEBUG: print(f"[I] {seq}")
        iterations += 1
    # check if length of seq is too long & homopolymer exists
    while len(seq) > desired_length:
        # INSERTION ERROR:
        # extract cumulative error position for all homopolymer locations
        possible_insertions = [res[i] if results[2][i] > 0.0 else 0.0 for i in range(len(res))]
        seq = seq[:np.argmax(possible_insertions)] + seq[np.argmax(possible_insertions) + 1:]
        if DEBUG: print(f"[D] {seq}")
        no_changes += 1  # easy since we always remove one base
        # re-calculate error values:
        results = calc_errors(seq)
        res = np.array([sum(x) for x in zip(*results)], dtype=np.float64)
    # if no homopolymer exists but sequence is too long, try heuristic repair but tag as low quality

    while sum([sum(x) for x in zip(*results[:-1])]) > 0 and iterations < 100:
        # MUTATION ERROR:
        # max_arg = np.random.choice(allmax(res))
        max_arg = allmax(res)
        max_arg = max_arg[int(len(max_arg) / 2)]
        possible_bases = allowed_bases - {seq[max_arg]}
        min_seq = seq
        min_sum_res = sum(res)
        for base in possible_bases:
            seq = seq[:max_arg] + base + seq[max_arg + 1:]
            results = calc_errors(seq)
            res = [sum(x) for x in zip(*results[:-1])]
            if sum(res) < min_sum_res:
                min_seq = seq
                min_sum_res = sum(res)
            if sum(res) == 0:
                break
        results[4][max_arg] = -min_sum_res  # we do not want to touch the same base again...
        res = [sum(x) for x in zip(*results)]
        no_changes += 1
        seq = min_seq
        if DEBUG: print(f"[M] {seq}")
        iterations += 1
    return seq, no_changes


def repair_all_n_chars(seq):
    # TODO: try to repair each "N" in the current sequence and set the quality to a low value (unless only one option is available)
    pass


def quality_score_to_phred_error_prob(quality_score, quality_score_format="Illumina 1.8+"):
    phred_offset = quality_score_format_dict[quality_score_format]
    quality_score = quality_score - phred_offset
    return 10 ** -(quality_score / 10) # accurecy = 1 - quality_score_to_phred_error_prob(quality_score)


# - try repair without brute force
#   -> create constraint combinations
# - use sequence length as a constraint
#   -> if longer than certain length & homopolymer: try removing bases from homopolymer

# TODO:
# - ideally we would want to check if there is a "close" sequence that adheres to all restrictions (clustering)
if __name__ == "__main__":
    import dinopy

    fqr_1 = dinopy.FastqReader("/home/michael/Code/RepairNatrix/mosla_logo/MOSLA1_A_R1.fastq.gz")
    print(fqr_1)
    print(fqr_1.guess_quality_format())
    for i, seq in enumerate(fqr_1.reads(True, True, True)):
        #print(f"{i}: {seq} : {min(seq.quality)}, {max(seq.quality)} {quality_score_to_phred_error_prob(max(seq.quality))}")
        print(f"{[quality_score_to_phred_error_prob(x) for x in seq.quality]}")

    DEBUG = False
    sequence = "ACCCAGTGACCCTATAGCCCCGTAGCTGCCACCCACCCACCCCCCCAAA"
    # sequence = "ACCCAGTGACCCTATAGCCCGTAGCTGCCACCCACCCACCCAAA"
    print(sequence)
    print(calc_errors(sequence))
    new_seq = try_repair(sequence, 44)
    print(new_seq)
    print(calc_errors(new_seq[0]))
