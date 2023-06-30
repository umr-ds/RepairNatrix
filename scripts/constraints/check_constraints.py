import functools
import glob
import os.path
from multiprocessing import freeze_support
from shutil import move
import json
import math
import typing
from collections import namedtuple

__FASTQ_EXT = ".fastq"

try:
    # we add these imorts to make the IDE happy about "snakemake.*" not existing...
    import failed_import, snakemake
except:
    pass

import logging
import multiprocessing
from tqdm import tqdm
from functools import partial
import numpy as np

pbar: tqdm = None

# TODO: use cdnarules where possible
import cdnarules
from FastDNARules import FastDNARules

rules = FastDNARules()
overall_gc_content_error_val = lambda data, low, high: [rules.overall_gc_content(data, calc_func=lambda x: 1.0 if (
        x < low or x > high) else 0.0)] * len(data)
windowed_gc_content_error_val = lambda data, window, low, high: [rules.windowed_gc_content(data, window_size=window,
                                                                                           ignore_last=True,
                                                                                           calc_func=lambda x: 1.0 if (
                                                                                                   x < low or x > high) else 0.0)] * len(
    data)

from homopolymers import homopolymer_error_val
from kmer import kmer_counting_error_val
from undesired_subsequences import UndesiredSubSequenceFinder

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


def quality_score_to_phred_error_prob(quality_score, quality_score_format="Illumina 1.8+"):
    phred_offset = quality_score_format_dict[quality_score_format]
    quality_score = quality_score - phred_offset
    return 10 ** -(quality_score / 10)


def phred_error_prob_to_quality_score(phred_error_prob, quality_score_format="Illumina 1.8+"):
    phred_offset = quality_score_format_dict[quality_score_format]
    return -10 * np.log10(phred_error_prob) + phred_offset


def read_fasta(filename: str) -> typing.Dict:
    """

    :param filename:
    :return: iterator of tuples: (name, seqeunce)
    """
    fasta_dict = {}
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('>'):
                i = 0
                seq_name = line.strip().split()[0][1:] + str(i)
                while seq_name in fasta_dict:
                    i += 1
                    seq_name = line.strip().split()[0][1:] + str(i)
                fasta_dict[seq_name] = ''
            else:
                fasta_dict[seq_name] += line.strip()
    return fasta_dict


def read_fastq(filename: str) -> typing.Iterator:
    """

    :param filename:
    :return: iterator of tuples: (name, sequence, phread quality score)
    """
    fastq_list = []
    last_line = -1
    with open(filename, 'r') as f:
        for j, line in enumerate(f):
            if line.startswith('@'):
                if last_line == 3:
                    fastq_list.append((seq_name, dna_seq, comment_line, phred_score))
                seq_name = line.strip()[1:]
                last_line = 0
            elif j % 4 == 1 and last_line == 0:
                dna_seq = line.strip()
                last_line = 1
            elif line.startswith('+') and last_line == 1:
                comment_line = line.strip()[1:]
                last_line = 2
            elif j % 4 == 3 and last_line == 2:
                phred_score = line.strip()
                last_line = 3
        if last_line == 3:
            fastq_list.append((seq_name, dna_seq, comment_line, phred_score))
    return iter(fastq_list)


def write_fasta(filename: str, data):
    """
    : typing.Union[
    typing.List[typing.Tuple[str, str]], typing.Iterable[typing.NamedTuple[typing.Union[str, bytes], str]]]
    :param filename:
    :param data: data in the format: List(Tuple(name, sequence))
    :return:
    """
    with open(filename, "w") as fp:
        for elem in data[:-1]:
            fp.write(f">{elem[0]}\n{elem[1]}\n")
        fp.write(f">{data[-1][0]}\n{data[-1][1]}")


def write_fastq(filename: str, data: typing.Union[
    typing.List[typing.Tuple[str, str, str]], typing.Tuple[typing.Union[str, bytes], str, str]]):
    """
    :param filename:
    :param data:
    :return:
    """
    with open(filename, "w") as fp:
        for elem in data[:-1]:
            fp.write(f">{elem[0]}\n{elem[1]}\n{elem[2]}\n{elem[3]}\n")
        fp.write(f">{elem[0]}\n{elem[1]}\n{elem[2]}\n{elem[3]}")


try:
    logging.basicConfig(filename=str(snakemake.log),
                        level=logging.DEBUG)
    inplace_repair = snakemake.params.inplace_repair
    repair_quality_score = snakemake.params.repair_quality_score
    allowed_min_gc = snakemake.params.min_gc_content
    allowed_max_gc = snakemake.params.max_gc_content
    allowed_max_gc_window_size = snakemake.params.windowed_gc_content_max_gc_window_size
    allowed_window_min_gc = snakemake.params.windowed_gc_content_min_gc_content
    allowed_window_max_gc = snakemake.params.windowed_gc_content_max_gc_content
    allowed_max_homopolymer_length = snakemake.params.max_homopolymer_length
    kmer_active = snakemake.params.kmer_active
    kmer_k = snakemake.params.kmer_k
    kmer_max_count = snakemake.params.kmer_max_count
    cores = snakemake.threads
    undesired_subsequences_file = snakemake.params.subsequences_file
    length = snakemake.params.sequence_length
    MAX_ITERATIONS = snakemake.params.maximum_repair_cycles
    USE_QUALITY_MAPPING = snakemake.params.use_quality_mapping
    snakemake_input_file_zero = snakemake.input[0]
except NameError as ne:
    logging.basicConfig(level=logging.DEBUG)
    logging.warning("No snakemake detected. Using default values for parameters.")
    inplace_repair = False
    repair_quality_score = 33
    allowed_min_gc = 0.3
    allowed_max_gc = 0.7
    allowed_max_gc_window_size = 50
    allowed_window_min_gc = 0.4
    allowed_window_max_gc = 0.6
    allowed_max_homopolymer_length = 4
    kmer_active = True
    kmer_k = 10
    kmer_max_count = 20
    cores = 8
    length = 156
    undesired_subsequences_file = "undesired_subsequences.txt"
    MAX_ITERATIONS = 50
    USE_QUALITY_MAPPING = True
    snakemake_input_file_zero = "../../results/assembly/40Dorn429_A/40Dorn429_A_cluster.fasta"

if allowed_min_gc < 1:
    allowed_min_gc = allowed_min_gc * 100
if allowed_max_gc < 1:
    allowed_max_gc = allowed_max_gc * 100
if allowed_window_min_gc < 1:
    allowed_window_min_gc = allowed_window_min_gc * 100
if allowed_window_max_gc < 1:
    allowed_window_max_gc = allowed_window_max_gc * 100
logging.info(
    f"Repairing file: {snakemake_input_file_zero} {'using in-place repair' if inplace_repair else 'adding repaired sequences.'}")
undesired_subsequences_finder = UndesiredSubSequenceFinder(undesired_subsequences_file)

# create quality mapping dict to speedup potential dereplication and repair cycles
quality_mapping_file = snakemake_input_file_zero.replace("_assembled.fastq", "derep.fastq")
quality_mapping_dict = dict()

if USE_QUALITY_MAPPING:
    if not os.path.exists(quality_mapping_file):
        # we do NOT have a dereplicated fastq file, thus we have to extract this info from the assembly file
        quality_mapping_file = snakemake_input_file_zero
        quality_mapping_reads = read_fastq(quality_mapping_file)
        for name, sequence, comment, quality in quality_mapping_reads:
            if sequence not in quality_mapping_dict.keys() \
                    or min(quality_mapping_dict[sequence]) < min(quality) \
                    or np.average(
                [int(quality_score_to_phred_error_prob(ord(x))) for x in quality_mapping_dict[sequence]]) \
                    < np.average([int(quality_score_to_phred_error_prob(ord(x))) for x in quality]):
                quality_mapping_dict[sequence] = quality
    else:
        quality_mapping_reads = read_fastq(quality_mapping_file)
        # file is already dereplicated
        for name, sequence, comment, quality in quality_mapping_reads:
            quality_mapping_dict[sequence] = quality


def repair_single_cluster(single_cluster_data, desired_length=160):
    global pbar
    pbar.update(1)
    centroid, cluster = single_cluster_data
    all_violations = calc_errors(centroid)
    res = [sum(x) for x in zip(*all_violations)]
    if sum(res) == 0 and len(res) == desired_length:
        # there was no error in the centroid. OPTIONAL: mark as correct
        with open("org_vs_repair.csv", "a") as fp:
            fp.write(f"{centroid}, {centroid}\n")
        return "C", centroid  # Centroid was Correct
        # continue
    else:
        # there was an error in the centroid: iterate over the cluster
        #   and find the sequence that fulfills all constraints and is closest to the centroid
        found = False
        for seq in cluster:
            seq_violations = calc_errors(seq)
            res = [sum(x) for x in zip(*seq_violations)]
            if sum(res) == 0 and len(res) == desired_length:
                with open("org_vs_repair.csv", "a") as fp:
                    fp.write(f"{centroid}, {seq}\n")
                return "S_C", seq  # Centroid was substituted - new centroid is correct - no repair was needed!
            # else:
            # the sequence does not fulfill all constraints
            # continue
        if not found:
            # no sequence in the cluster fulfills all constraints
            #   -> repair the centroid
            possible_results = []
            pair_with_quality = namedtuple("cluster_repair", ["name", "seq", "quality"])
            for i, elem in enumerate([centroid] + cluster):
                q_score = quality_mapping_dict.get(elem.upper(), repair_quality_score)
                repair_res = try_repair(pair_with_quality(f"cluster_repair_{i}", elem, q_score), desired_length)
                repair_violations = calc_errors(repair_res[1])
                res = [sum(x) for x in zip(*repair_violations)]
                if sum(res) == 0 and len(res) == desired_length:
                    # run repair for all elements in the cluster
                    # and take the repaired sequence with the LEAST repairs needed!
                    # centroid was substituted - substituted sequence had to be repaired
                    possible_results.append((repair_res[4], (f"S_R_C_{repair_res[4]}", repair_res[1])))
                else:
                    continue
            if len(possible_results) > 0:
                with open("org_vs_repair.csv", "a") as fp:
                    fp.write(f"{centroid}, {sorted(possible_results, key=lambda x: x[0])[0][1][1]}\n")
                return sorted(possible_results, key=lambda x: x[0])[0][1]
        with open("org_vs_repair.csv", "a") as fp:
            fp.write(f"{centroid}, _\n")
        return "F", centroid


def repair_clusters(desired_length):
    """
    clusters: list of Clusters, each cluster has a centroid and a list of sequences
    (with a distance-value to the centroid)
    IMPORTANT: clusters are sorted by distance to centroid!
    desired_length: int, desired length of the sequences in the clusters
    :returns: list of centroids with a quality score indicating the performed repairs together with the initial quality
    """
    global pbar
    try:
        input_file = snakemake.input.cent  # centroid file...
    except NameError:
        logging.info("Running in non-snakemake mode - this should be used for testing only")
        input_file = snakemake_input_file_zero
    logging.info(f"Parsing cluster for input file: {input_file}")

    clusters = []

    # get folder for a file string:
    input_folder = os.path.dirname(input_file)
    for cluster_file in glob.glob(f"{input_folder}/clust*"):
        cluster = [seq for name, seq in read_fasta(cluster_file).items()]
        clusters.append((cluster[0], cluster[1:]))
    clusters = sorted(clusters, key=lambda x: len(x[1]), reverse=True)
    output_file = input_file.replace("cluster.fasta", "repaired.fasta")
    output_cluster_mapping_file = input_file.replace("cluster.fasta", "repaired_clusters.json")
    # run different repairs based on the input-file:
    # if we have a list of clusters(+centroids) we want to run a special repair for each cluster:
    # - if the centroid does fulfill all constraints: return centroid as correct
    # - if the centroid does not fulfill all constraints:
    #     - search in the cluster for the sequence that is closest to the centroid and fulfills all constraints
    #     - if no such sequence exists:
    #         - first: repair only the centroid
    #         - if this fails after CHANGE_LIMIT changes:
    #           - repair _ALL_ sequences in the cluster and choose the one with the least changes
    #             _AND_ the shortest distance to the centroid
    # res_centroids = []
    pbar = tqdm(total=len(clusters))
    p = multiprocessing.Pool(cores)
    # res_centroids = [x for x in
    #                 p.imap_unordered(partial(repair_single_cluster, desired_length=desired_length), iterable=clusters,
    #                                  chunksize=max(1,math.ceil(len(clusters) / (cores * 10))))]
    res_centroids = [partial(repair_single_cluster, desired_length=desired_length)(y) for y in clusters]
    if not inplace_repair:
        # read all entries of input_file (containing all initial centroids)
        # initial_centroids = dinopy.FastqReader(input_file)
        res_seqs = set([b for a, b in res_centroids])
        # add original centroid to
        print(f"{input_file}")
        for _name, _sequence in read_fasta(input_file).items():
            if _sequence not in res_seqs:
                res_centroids.append((f"{_name.replace('_', '-')}_O_F", _sequence))

    if os.path.exists(output_file):
        renamed_file = output_file.replace(__FASTQ_EXT, "old_.fastq").replace(".fasta", "_old.fasta")
        logging.warning(f"[WARNING] File already exists, renaming old file to: {renamed_file}")
        move(output_file, renamed_file)
    write_fasta(output_file, [(str(i) + "_" + a, b) for i, (a, b) in
                              enumerate(sorted(res_centroids, key=lambda tpl: sort_results(tpl[0])))])
    if os.path.exists(output_cluster_mapping_file):
        renamed_file = output_cluster_mapping_file.replace(".json", "_old.json")
        logging.warning(f"[WARNING] File already exists, renaming old file to: {renamed_file}")
        move(output_cluster_mapping_file, renamed_file)
    with open(output_cluster_mapping_file, "w") as cluster_output:
        json.dump([x for x in zip([b for a, b in res_centroids], clusters)], fp=cluster_output)
    return res_centroids


def sort_results(in_str):
    if in_str is None:
        return MAX_ITERATIONS + 5000
    try:
        # ensure we have str and not bytes
        in_str = in_str.decode()
    except:
        pass
    if ":" in in_str:
        in_str = in_str.split(":")[-1]
        in_str = in_str[in_str.find("_") + 1:]
    if in_str == "C":  # centroid was correct
        return -1000
    elif in_str == "S_C":  # centroid was substituted from cluster, new centroid correct
        return -500
    elif in_str.startswith("S_R_C_"):  # centroid was (potentially) substituted; new centroid was repaired (X repairs)
        return int(in_str.split("_")[-1])
    elif in_str == "F":  # repair failed (repair limit reached for all elements in the cluster!)
        return MAX_ITERATIONS + 1000
    elif in_str.startswith("O_F"):  # in case of INPLACE_REPAIR = False, we want to flag invalid original seqs as failed
        return MAX_ITERATIONS + 1001
    else:
        logging.error(f"sorting found not handled case: {in_str}")
        return MAX_ITERATIONS


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
               windowed_gc_content_error_val(seq, int(allowed_max_gc_window_size), float(allowed_window_min_gc),
                                             float(allowed_window_max_gc)),
               homopolymer_error_val(seq, int(allowed_max_homopolymer_length), True),
               kmer_counting_error_val(seq, int(kmer_k), int(kmer_max_count), kmer_active),
               [0.0] * len(seq)]
    return results


def try_repair(seq_instance: typing.NamedTuple, desired_length, update_pbar=False):
    if update_pbar:
        pbar.update(1)
    if "quality" not in seq_instance or seq_instance[2] is None:
        # we might be in a FASTA-File and thus have no access to the quality scores
        phred_score = [repair_quality_score] * len(seq_instance[1])
    else:
        phred_score = seq_instance[2]

    if not isinstance(seq_instance[1], str):
        seq = seq_instance[1].decode("utf-8")
    else:
        seq = seq_instance[1]
    org_seq = seq
    iterations = 0
    no_changes = 0
    results = calc_errors(seq)
    res = [sum(x) for x in zip(*results)]
    assert len(phred_score) == len(seq)
    check_and_modified_positions = set()
    while len(seq) < desired_length and no_changes < MAX_ITERATIONS:
        # DELETION ERROR:
        # this is a bit more complicated: we cant simply pick the position with the highest error-value
        # because errors might appear _after_ a possible deletion
        # however we _can_ be sure that the deletion should be _before_

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
            # if the error value does not change, we should take the base that would more likely repair in the long run:
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
            no_changes = 0
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
            # if the gc-content is off:
            # choose the _first_ subsequence of correct length that has a correct gc-content
            # OR: use a sliding window to find the best subsequence:
            # first or last window should be choosen (either there was an offset at the beginning or
            # there were additional bases at the end)

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
    return org_seq, seq, phred_score, seq_instance[0], no_changes


def main(desired_length=160):
    global pbar
    repaired_tuples = []
    phred_error_prob_to_quality_score(quality_score_to_phred_error_prob(30))
    # Read in the fasta file _AFTER_ clustering
    try:
        input_file = snakemake.input[0]
    except NameError:
        logging.info("Running in non-snakemake mode - this should be used for testing only")
        input_file = snakemake_input_file_zero
    logging.info(f"Input files: {input_file}")

    if input_file.endswith(__FASTQ_EXT):
        seqs: typing.List[typing.List] = [[name, sequence, quality] for name, sequence, comment, quality in
                                          read_fastq(input_file)]
    else:
        logging.info("Fasta-File - No phred score given. Using default value.")
        seqs: typing.List[typing.List] = [[name, sequence, repair_quality_score] for name, sequence in
                                          read_fasta(input_file).items()]
    pbar = tqdm(total=len(seqs))
    p = multiprocessing.Pool(cores)
    a = [x for x in p.imap_unordered(partial(try_repair, desired_length=desired_length, update_pbar=True),
                                     iterable=seqs,
                                     chunksize=math.ceil(len(seqs) / (cores * 10)))]
    out_data = []
    if not inplace_repair:
        for seq in seqs:
            seq_name = seq[0].decode() + "_org"
            if "quality" not in seq or seq[2] is None:
                # we might be in a FASTA-File and thus have no access to the quality scores
                phred_score = [int(repair_quality_score)] * len(seq[1])
            else:
                phred_score = [int(x) for x in seq[2]]
            out_data.append((seq[1], seq_name.encode(), [chr(int(x)) for x in phred_score]))
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
        except AttributeError:  # support for non-snakemake runs
            res = elem[1]
        out_data.append((res.encode("utf-8"), seq_name.encode(), bytes_quality_values))
    out_path = input_file.replace(__FASTQ_EXT, "_constraint_repaired.fastq").replace(".fasta",
                                                                                     "_constraint_repaired.fastq")
    if os.path.exists(out_path):
        renamed_file = out_path.replace(__FASTQ_EXT, "_old.fastq").replace(".fasta", "_old.fasta")
        logging.warning(f"[WARNING] File already exists, renaming old file to: {renamed_file}")
        move(out_path, renamed_file)
    write_fastq(out_path, out_data)
    logging.info(f"Saving result: {out_path}")

    out_path = input_file.replace(__FASTQ_EXT, "_constraint_repaired_mapping.json").replace(".fasta",
                                                                                            "_constraint_repaired_mapping.json")
    if os.path.exists(out_path):
        renamed_file = out_path.replace(".json", "_old.json")
        logging.warning(f"[WARNING] File already exists, renaming old file to: {renamed_file}")
        move(out_path, renamed_file)
    with open(out_path, "w") as f:
        json.dump(repaired_tuples, f)
    logging.info(f"Saving result-mapping for alternative centroid selection during decoding: {out_path}")


if __name__ == '__main__':
    freeze_support()
    repair_clusters(length)
