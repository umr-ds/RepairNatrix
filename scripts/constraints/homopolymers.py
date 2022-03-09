import random
import typing

try:
    import cdnarules
except ModuleNotFoundError:
    print("C Module failed to load, falling back to slow mode")


def homopolymer_old(sequence, count):
    """ edge-case missing: if homopolymer is at the end of the sequence! """
    prev_char = sequence[0]
    # we add a symbol at the end to catch edge-case with the homopolymer beeing at the end of the sequence:
    sequence = sequence + "_"
    curr_start_pos = 0
    length = len(sequence)
    for char_pos in range(1, length):
        if prev_char != sequence[char_pos]:
            homopolymer_length = char_pos - curr_start_pos
            if homopolymer_length >= count:
                return True
            curr_start_pos = char_pos
            prev_char = sequence[char_pos]
    return False


def homopolymer_val(sequence):
    res = [0.0] * len(sequence)
    prev_char = sequence[0]
    # we add a symbol at the end to catch edge-case with the homopolymer beeing at the end of the sequence:
    sequence = sequence + "_"
    curr_start_pos = 0
    length = len(sequence)
    for char_pos in range(1, length):
        if prev_char != sequence[char_pos]:
            homopolymer_length = char_pos - curr_start_pos
            res[curr_start_pos:char_pos] = [1.0 * homopolymer_length] * (char_pos - curr_start_pos)
            curr_start_pos = char_pos
            prev_char = sequence[char_pos]
    return res


def homopolymer_error_val(sequence, count, increasing_chance=False):
    return [(x-count if increasing_chance else 1.0) if x > count else 0.0 for x in homopolymer_val(sequence)]


def homopolymer(data, count):
    """
    Calculates the dropchance based on the chance of homopolymers to mutate. Homopolymers are repeats of the same
    single unit and have a higher chance to mutate than regular polymers.
    :param data: The DNA sequence to check for homopolymers.
    :return: The dropchance based on the occurence of homopolymers.
    """
    return longestSequenceOfChar_python(data, '*')[1] >= count


def longestSequenceOfChar_python(text: typing.AnyStr, char_x="*") -> typing.Tuple[str, int]:
    n = len(text)
    c = 0
    res = char_x
    curr = 1
    i = 0
    while i < n:
        if i < n - 1 and text[i] == text[i + 1]:
            curr += 1
        else:
            if curr > c and (text[i] == char_x or char_x == "*"):
                c = curr
                res = text[i]
            curr = 1
        i += 1
    return res, c


def longestSequenceOfChar(text: typing.AnyStr, char_x="*") -> typing.Tuple[str, int]:
    try:
        return cdnarules.longestSequenceOfChar(text, char_x)
    except NameError:
        print("error...")
        return longestSequenceOfChar_python(text, char_x)


if __name__ == "__main__":
    print(homopolymer("AAAAAAAAAACACACTTTTTTTTAAAAAAAAAAA", 5))
    print(homopolymer("AAAAGGGGGGCGGGGAGCCCTTTTTCGCGCCCCCCGGGGTTTTTT", 5))
    print(homopolymer("AACCGGACAGGGGG", 5))
    print(not homopolymer("AACCGGACAGGGG", 5))
    print(homopolymer_val("AAAAAAAAAACACACTTTTTTTTAAAAAAAAAAA"))
    print(homopolymer_error_val("AAAAAAAAAACACACTTTTTTTTAAAAAAAAAAACCCCC", 5, True))
    for _ in range(200):
        a = "".join([random.choice(["A", "C", "G", "T"]) for _ in range(100)])
        # print(a)
        assert (homopolymer(a, 5) == homopolymer_old(a, 5))
    import timeit

    a = "".join([random.choice(["A", "C", "G", "T"]) for _ in range(100)])
    res_a = timeit.timeit("""a = "".join([random.choice(["A", "C", "G", "T"]) for _ in range(100)])
homopolymer(a, 5)""", setup="from __main__ import homopolymer,homopolymer_old, random, cdnarules",
                          number=10000)
    #print(f"new: {res_a=}")
    res_b = timeit.timeit("""a = "".join([random.choice(["A", "C", "G", "T"]) for _ in range(100)])
homopolymer_old(a, 5)""",
                          setup="from __main__ import homopolymer_old,homopolymer, random, cdnarules",
                          number=10000)
    #print(f"old: {res_b=}")
