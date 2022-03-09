import regex as re


def undesired_subsequences_val(sequence, file):
    res = [0.0] * len(sequence)
    with open(file) as f:
        list_of_subsequences = [seq.strip() for seq in f.readlines()]

    regextext = "|".join(list_of_subsequences)
    for m in re.finditer(regextext, sequence, overlapped=True):
        res[m.span()[0]:m.span()[1]] = [1.0] * (m.span()[1] - m.span()[0])
    return res

def undesired_subsequences(sequence, file):
    return max(undesired_subsequences_val(sequence, file)) > 0.0

if __name__ == "__main__":
    print(undesired_subsequences("ACCCAGTGACCCTATAGCCCCGTAGCTGCCACCCACCCACCCCAAA", "undesired_subsequences.txt"))
    print(undesired_subsequences_val("ACCCAGTGACCCTATAGCCCCGTAGCTGCCACCCACCCACCCCAAA", "undesired_subsequences.txt"))