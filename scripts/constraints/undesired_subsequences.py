import regex as re

def undesired_subsequences(sequence, file):
    with open(file) as f:
        list_of_subsequences = [seq.strip() for seq in f.readlines()]

    regextext = "|".join(list_of_subsequences)
    for m in re.finditer(regextext, sequence, overlapped=True):
        return True
    return False