from math import ceil

def windowed_gc_content(sequence, gc_min, gc_max, window_size):
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


def overall_gc_content(sequence, gc_min, gc_max):
    length = len(sequence)
    basecount = dict()
    for char_pos in range(length):
        if sequence[char_pos] in basecount:
            basecount[sequence[char_pos]] += 1
        else:
            basecount[sequence[char_pos]] = 1
    gc_sum = 1.0 * (basecount["G"] if "G" in basecount else 0.0 + basecount["C"] if "C" in basecount else 0.0)
    gc_percentage = 1.0 * (gc_sum / length * 100) / 100.0
    if not (gc_min <= gc_percentage <= gc_max):
        return True
    return False


if __name__ == "__main__":
    print(overall_gc_content("AAAAAAAAGAACACACTTTTTTTTAAAAAAAAAAA",0.3,0.7))
    print(windowed_gc_content("GCGCATATATATATATATAGCGCGCGCGCGCGCG",0.3,0.7,10))
