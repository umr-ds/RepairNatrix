def default_error_function(homopolymer_length):
    return 1.0 if homopolymer_length > 5 else 0.0


def __old_default_error_function(homopolymer_length, base=None):
    """
    Default function to calculate error probabilities based on the length of homopolymers.
    :param homopolymer_length:
    :param base:
    :return: Error probability of the homopolymer.
    """
    if homopolymer_length < 3:
        return 0
    elif homopolymer_length < 4:
        return 0.3
    elif homopolymer_length < 5:
        return 0.6
    elif homopolymer_length < 6:
        return 0.9
    else:
        return 1.0


def create_result(startpos: int, endpos: int, base: chr):
    """
    Creates a dictionary containing all results.
    :param startpos:
    :param endpos:
    :param errorprob:
    :param i:
    :return:
    """
    res = {'start': startpos, 'end': endpos, 'base': base}
    return res


def homopolymer(sequence, apply_for=None, error_function=default_error_function):
    """
    Generates a list of dictionarys containing base, startpos, endpos and errorprobability for all homopolymers in
    the given sequence.
    :param sequence:
    :param apply_for:
    :param error_function:
    :return:
    """
    if error_function is None:
        error_function = default_error_function
    result = []
    if apply_for is None:
        apply_for = {'G', 'T', 'A', 'C'}
    prev_char = sequence[0]
    curr_start_pos = 0
    length = len(sequence)
    for char_pos in range(1, length):
        if prev_char != sequence[char_pos]:
            error_prob = error_function(char_pos - curr_start_pos)
            if error_prob > 0.0:
                result.append(create_result(curr_start_pos, char_pos - 1, prev_char))
            curr_start_pos = char_pos
            prev_char = sequence[char_pos]
    error_prob = error_function(length - curr_start_pos)
    if error_prob > 0.0:
        result.append(create_result(curr_start_pos, length - 1, prev_char))
    return result


if __name__ == "__main__":
    print(homopolymer("AAAAAAAAAACACACTTTTTTTTAAAAAAAAAAA"))
    print(homopolymer("AAAAGGGGGGCGGGGAGCCCTTTTTCGCGCCCCCCGGGGTTTTTT"))
