from math import ceil
from collections import Counter


def windowed_gc_content_old(sequence, gc_min, gc_max, window_size):
    length = len(sequence)
    no_windows = ceil(1.0 * length / window_size)
    for i in range(no_windows):
        basecount = dict()
        window_sequence = sequence[i * window_size:min((i + 1) * window_size, length)]
        curr_length = len(window_sequence)
        for char_pos in range(curr_length):
            if window_sequence[char_pos] in basecount:
                basecount[window_sequence[char_pos]] += 1
            else:
                basecount[window_sequence[char_pos]] = 1
        gc_sum = 1.0 * ((basecount["G"] if "G" in basecount else 0.0) + (basecount["C"] if "C" in basecount else 0.0))
        gc_percentage = 1.0 * (gc_sum / curr_length * 100.0) / 100.0
        if not (gc_min <= gc_percentage <= gc_max):
            return True
    return False


def windowed_gc_content(sequence, gc_min, gc_max, window_size):
    for i in range(0, len(sequence), window_size):
        if overall_gc_content(sequence[i:i + window_size], gc_min, gc_max, ):
            return True
    return False


def windows_gc_content_error_val(sequence, window_size, gc_min, gc_max):
    return [1.0 if (x < gc_min or x > gc_max) else 0.0 for x in windowed_gc_content_val(sequence, window_size)]


def windowed_gc_content_val(sequence, window_size):
    res = []
    for i in range(0, len(sequence), window_size):
        window_val = overall_gc_content_val(sequence[i:i + window_size])
        for _ in range(window_size):
            res.append(window_val)
    return res


def overall_gc_content(sequence, gc_min, gc_max):
    gc_percentage = overall_gc_content_val(sequence)
    if not (gc_min <= gc_percentage <= gc_max):
        return True
    return False


def overall_gc_content_val(sequence):
    counter = Counter(sequence)
    count = counter["G"] + counter["C"]
    return (count / len(sequence) * 100) / 100


def overall_gc_content_error_val(sequence, gc_min, gc_max):
    x = overall_gc_content_val(sequence)
    tmp = 1.0 if (x < gc_min or x > gc_max) else 0.0
    return [tmp for _ in range(len(sequence))]


if __name__ == "__main__":
    import random

    for _ in range(2):
        a = "".join([random.choice(["A", "C", "G", "T"]) for _ in range(100)])
        # print(a)
        print(windowed_gc_content(a, 0.3, 0.7, 10) == windowed_gc_content_old(a, 0.3, 0.7, 10))
    seq = "ACCCACACTGGGGGCCGCGCGCGCGCCGATATATATATCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGACTGTGCTAGTAGTAGATGATACC"
    print(windowed_gc_content_val(seq, 10))
    print(windows_gc_content_error_val(seq, 10, 0.3, 0.7))
    print(overall_gc_content_error_val(seq, 0.3, 0.7))