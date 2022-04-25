import os.path
from shutil import move
import json
import math
import typing
from collections import namedtuple

import dinopy
import logging
import multiprocessing
from functools import partial
import numpy as np

# try:
# TODO. use cdnarules where possible
import cdnarules
# except ImportError:
from gc_content import overall_gc_content_error_val
from homopolymers import homopolymer_error_val
from kmer import kmer_counting_error_val
from undesired_subsequences import UndesiredSubSequenceFinder

MAX_ITERATIONS = 50
allowed_bases = {"A", "C", "G", "T"}
DEBUG = False
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

try:
    logging.basicConfig(filename=str(snakemake.log),
                        level=logging.DEBUG)
    inplace_repair = snakemake.params.inplace_repair
    repair_quality_score = snakemake.params.repair_quality_score
    allowed_min_gc = snakemake.params.min_gc_content
    allowed_max_gc = snakemake.params.max_gc_content
    allowed_max_homopolymer_length = snakemake.params.max_homopolymer_length
    kmer_k = snakemake.params.kmer_k
    kmer_max_count = snakemake.params.kmer_max_count
    cores = snakemake.threads
    undesired_subsequences_file = snakemake.params.subsequences_file
except NameError as ne:
    logging.basicConfig(level=logging.DEBUG)
    logging.warning("No snakemake detected. Using default values for parameters.")
    inplace_repair = False
    repair_quality_score = 33
    allowed_min_gc = 0.4
    allowed_max_gc = 0.7
    allowed_max_homopolymer_length = 3
    kmer_k = 10
    kmer_max_count = 20
    cores = 8
    undesired_subsequences_file = "undesired_subsequences.txt"

undesired_subsequences_finder = UndesiredSubSequenceFinder(undesired_subsequences_file)

""" 
idea:
    - AFTER clustering: check if there are multiple clusters with the same seed (i.e. the same X starting bases)
        - parameter: unique first X bases (int)
    - BEFORE forward/reverse-read-merging: limit GC content rule to windowed mode and skip last/first window.
"""


def allmax(a, return_index=True):
    if len(a) == 0:
        return []
    all_ = [0]
    max_ = a[0]
    all_max_ = [a[0]]
    for i in range(1, len(a)):
        if a[i] > max_:
            all_ = [i]
            max_ = a[i]
            all_max_ = [a[i]]
        elif a[i] == max_:
            all_.append(i)
            all_max_.append(a[i])
    if return_index:
        return all_
    else:
        return all_max_


def calc_errors(seq):
    """
    Calculates the error-values for each constraint.
    :param seq: sequence to test
    :return: list of lists of error-values for each constraint:
        [undesired_subsequences, overall_gc_content, homopolymer_error, kmer_counting, list of 0 for minus optimization]
    """
    results = [undesired_subsequences_finder.undesired_subsequences_val(seq),
               overall_gc_content_error_val(seq, float(allowed_min_gc), float(allowed_max_gc)),
               homopolymer_error_val(seq, int(allowed_max_homopolymer_length), True),
               kmer_counting_error_val(seq, int(kmer_k), int(kmer_max_count)),
               [0.0] * len(seq)]
    return results


def repairInsertionWithDeletion(seq, desired_length):
    # the sequence is a result of (multiple) insertions and deletions
    # while this should not be the first option to choose (Occam's razor), this CAN happen...

    # first we enforce the desired length just like try_repair() does

    # then we try removing a base to lower the error value
    # -> if there is a different location with an error remaining, we then try to insert a base at that location
    # -> if there is no error remaining:
    #   -> the wrong base is removed
    #   -> the missing insertion did not create a rule-violation
    pass


def try_repair(seq_instance: typing.NamedTuple, desired_length):
    if "quality" not in seq_instance or seq_instance[2] is None:
        # we might be in a FASTA-File and thus have no access to the quality scores
        phred_score = [repair_quality_score] * len(seq_instance[1])
    else:
        phred_score = seq_instance[2]
    seq = seq_instance[1].decode("utf-8")
    org_seq = seq
    iterations = 0
    no_changes = 0
    results = calc_errors(seq)
    res = [sum(x) for x in zip(*results)]
    assert len(phred_score) == len(seq)
    check_and_modified_positions = set()
    while len(seq) < desired_length:
        # DELETION ERROR:
        # this is a bit more complicated: we cant simply pick the position with the highest error-value
        # because errors might appear _after_ a possible deletion
        # however we _can_ be sure that the deletion should be _before_

        # to/do: right now we randomly take the argmax, but we might want to go for the middle-position (minimizes the number of changes for homopolymers)
        # max_arg = np.random.choice(allmax(res))

        # find the positions with the highest error-value
        max_error_pos = allmax(res)
        if len(max_error_pos) == 0 or res[max_error_pos[0]] == 0:
            max_error_pos = [len(seq) - 1]  # we want to append at the end of the sequence...
        # from this set search for the positions with the lowest phred score
        max_arg = [max_error_pos[x] for x in allmax([-phred_score[y] for y in max_error_pos])]
        # filter out all positions that were already checked
        max_arg = [x for x in max_arg if x not in check_and_modified_positions]
        if len(max_arg) == 0:
            break
        # choose the middle position from the remaining positions
        max_arg = max_arg[int(len(max_arg) / 2)]

        possible_bases = allowed_bases - {seq[max_arg]}
        min_seq = seq
        min_phred_score = phred_score
        min_sum_res = sum(res)
        min_results = results
        min_res = res
        insert = 0
        for base in possible_bases:
            seq = seq[:max_arg] + base + seq[max_arg + insert:]  # for the first run we insert instead of replacing..
            phred_score = phred_score[:max_arg] + [repair_quality_score] + phred_score[max_arg + insert:]
            insert = 1
            results = calc_errors(seq)
            res = [sum(x) for x in zip(*results)]
            # TODO we might want to try out each base and choose the minimum error...
            # TODO: if the error value does not change, we should take the base that would more likely repair in the long run:
            #   -> we should select the base that moves the GC content closer to the allowed range
            if sum(res) < min_sum_res or len(check_and_modified_positions) == len(seq):
                min_seq = seq
                min_phred_score = phred_score
                min_sum_res = sum(res)
                min_results = results
                min_res = res
            if sum(res) == 0:
                break
        check_and_modified_positions.add(max_arg)
        if seq != min_seq:
            no_changes += 1
            seq = min_seq
            phred_score = min_phred_score
            results = min_results
            res = min_res
            # we need to shift all already checked positions that are bigger than the current base by one
            tmp = set()
            for pos in check_and_modified_positions:
                if pos <= max_arg:
                    tmp.add(pos)
                else:
                    tmp.add(pos + 1)
            check_and_modified_positions = tmp
            if DEBUG: print(f"[I] {seq}")
        iterations += 1
    # check if length of seq is too long & homopolymer exists
    while len(seq) > desired_length:
        # INSERTION ERROR:
        # extract cumulative error position for all homopolymer locations
        possible_insertions = [res[i] if results[2][i] > 0.0 else -1.0 for i in range(len(res))]
        if np.max(possible_insertions) == -1:
            #TODO: if the gc-content is off:
            # choose the _first_ subsequence of correct length that has a correct gc-content
            # OR: use a sliding window to find the best subsequence:
            # first or last window should be choosen (either there was an offset at the beginning or
            # there were additional bases at the end)

            # no homopolymer found, we can only remove the last base

            seq = seq[:-1]
            phred_score = phred_score[:-1]
        else:
            # since we are only removing homopolymers we can take ANY position with the highest error
            seq = seq[:np.argmax(possible_insertions)] + seq[np.argmax(possible_insertions) + 1:]
            phred_score = phred_score[:np.argmax(possible_insertions)] + phred_score[
                                                                         np.argmax(possible_insertions) + 1:]
        if DEBUG: print(f"[D] {seq}")
        no_changes += 1  # easy since we always remove one base
        # re-calculate error values:
        results = calc_errors(seq)
        res = np.array([sum(x) for x in zip(*results)], dtype=np.float64)
        assert len(phred_score) == len(seq)
    # if no homopolymer exists but sequence is too long, try heuristic repair but tag as low quality

    while sum([sum(x) for x in zip(*results[:-1])]) > 0 and iterations < MAX_ITERATIONS:
        # MUTATION ERROR:
        # max_arg = np.random.choice(allmax(res))
        # take all positions with the MOST errors, from this set take all sequences with the lowest phred score
        max_arg = allmax([-phred_score[x] for x in allmax(res)])
        max_arg = max_arg[int(len(max_arg) / 2)]
        possible_bases = allowed_bases - {seq[max_arg]}
        min_seq = seq
        min_sum_res = sum(res)
        for base in possible_bases:
            seq = seq[:max_arg] + base + seq[max_arg + 1:]
            phred_score = phred_score[:max_arg] + [repair_quality_score] + phred_score[max_arg + 1:]
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
    return (org_seq, seq, phred_score, seq_instance[0], no_changes)


def repair_all_n_chars(seq):
    # TODO: try to repair each "N" in the current sequence and set the quality to a low value (unless only one option is available)
    pass


def quality_score_to_phred_error_prob(quality_score, quality_score_format="Illumina 1.8+"):
    phred_offset = quality_score_format_dict[quality_score_format]
    quality_score = quality_score - phred_offset
    return 10 ** -(quality_score / 10)  # accurecy = 1 - quality_score_to_phred_error_prob(quality_score)


def phred_error_prob_to_quality_score(phred_error_prob, quality_score_format="Illumina 1.8+"):
    phred_offset = quality_score_format_dict[quality_score_format]
    return -10 * np.log10(phred_error_prob) + phred_offset


# - try repair without brute force
#   -> create constraint combinations
# - use sequence length as a constraint
#   -> if longer than certain length & homopolymer: try removing bases from homopolymer

# TODO:
# - ideally we would want to check if there is a "close" sequence that adheres to all restrictions (clustering)
def main():
    repaired_tuples = []
    phred_error_prob_to_quality_score(quality_score_to_phred_error_prob(30))
    # Read in the fasta file _AFTER_ clustering
    try:
        input_file = snakemake.input[0]
    except NameError:
        logging.info("Running in non-snakemake mode - this should be used for testing only")
        input_file = "../../results/assembly/MOSLA12_A/MOSLA12_A_assembled.fastq"
        # input_file = "/home/michael/Code/RepairNatrix/results/assembly/MOSLA2_A/MOSLA2_A_cluster.fasta"
    logging.info(f"Input files: {input_file}")
    # fasta_1 = dinopy.FastaReader(input_file)
    # print(try_repair(
    #    "GGATTGAGCAGCTGTCTATAGTACGTACTTTCAGATATTGCGATAAGCGTGTTGACCAAAACTTGCTGACGCATCGAGGAAAGATGTACTTCTTGGGGTCAGAGGGCTACCCTGATAGTTTTACGCAAGATCTCTCCGACGTCTGGTAATATTTTGCCGATGCAAAGTGATCTCAAGGCAGCGCACCCGCACTGCTGGTACGCCGATCAGAATGAGACTGACCGATACAGTGTTAAGCGAGTGAC",
    #    None, 160))
    # a = [try_repair(x.sequence, None, 160) for x in fasta_1.reads(True, str)]

    # with multiprocessing.Pool(cores) as p:
    # seqs = [x.sequence for x in fasta_1.reads(True, dtype=str)]
    p = multiprocessing.Pool(cores)
    # a = [x for x in p.imap_unordered(partial(try_repair, phred_score=None, desired_length=160),
    #                                 iterable=seqs,
    #                                 chunksize=math.floor(len(seqs) / (cores * 2)))]
    # print(a)
    # a = [x for x in map(partial(try_repair, phred_score=None, desired_length=160), seqs)]

    # possible repair: after assembly
    """
    # fqr_1 = dinopy.FastaReader("/home/michael/Code/RepairNatrix/results/assembly/MOSLA2_A/MOSLA2_A_derep.fasta")
    fqr_1 = dinopy.FastqReader(input_file)
    seqs = [x for x in fqr_1.reads()]
    try_repair(seqs[0], 160)
    """
    # better: possible repair: after dereplication


    if input_file.endswith(".fastq"):
        fqr_1 = dinopy.FastqReader(input_file)
        seqs: typing.List[typing.List] = [[read.name, read.sequence, read.quality] for read in fqr_1.reads()]
    else:
        logging.info("Fasta-File - No phred score given. Using default value.")
        fqr_1 = dinopy.FastaReader(input_file)
        seqs: typing.List[typing.List] = [[read.name, read.sequence, repair_quality_score] for read in
                                          fqr_1.reads(True)]

    a = [x for x in p.imap_unordered(partial(try_repair, desired_length=160),
                                     iterable=seqs,
                                     chunksize=math.ceil(len(seqs) / (cores * 10)))]
    # a = [x for x in map(partial(try_repair, desired_length=160), seqs)]
    out_data = []
    if not inplace_repair:
        for seq in seqs:
            seq_name = seq[0].decode() + "_org"
            if "quality" not in seq or seq[2] is None:
                # we might be in a FASTA-File and thus have no access to the quality scores
                phred_score = [int(phred_error_prob_to_quality_score(quality_score_to_phred_error_prob(repair_quality_score)))] * len(seq[1])
            else:
                phred_score = [int(phred_error_prob_to_quality_score(quality_score_to_phred_error_prob(x))) for x in seq[2]]
            out_data.append((seq[1], seq_name.encode("utf-8"), bytes(phred_score)))
    for elem in sorted(a, key=lambda x: x[-1]):
        # store result as a fastq file to reflect
        bytes_quality_values = bytes(
            [int(phred_error_prob_to_quality_score(quality_score_to_phred_error_prob(x))) for x in elem[2]])
        seq_name = elem[3].decode() + "_" + str(elem[4])
        # only add if the sequence was repaired
        if elem[4] > 0:
            repaired_tuples.append([elem[0], elem[1]])
        try:
            res = elem[1].decode("utf-8")
        except AttributeError: # support for non-snakemake runs
            res = elem[1]
        out_data.append((res.encode("utf-8"), seq_name.encode("utf-8"), bytes_quality_values))
    out_path = input_file.replace(".fastq", "_constraint_repaired.fastq").replace(".fasta", "_constraint_repaired.fastq")
    if os.path.exists(out_path):
        renamed_file = out_path.replace(".fastq", "_old.fastq")
        print(f"[WARNING] File already exists, renaming old file to: {renamed_file}")
        move(out_path, renamed_file)
    with dinopy.FastqWriter(out_path) as of:
        of.write_reads(out_data) # , dtype=str

    out_path = input_file.replace(".fastq", "_constraint_repaired_mapping.json").replace(".fasta", "_constraint_repaired_mapping.json")
    if os.path.exists(out_path):
        renamed_file = out_path.replace(".json", "_old.json")
        print(f"[WARNING] File already exists, renaming old file to: {renamed_file}")
        move(out_path, renamed_file)
    with open(out_path, "w") as f:
        json.dump(repaired_tuples, f)
    print(len([x for x in a if x[4] < 100]), len(a))


if __name__ == "__main__":
    main()
