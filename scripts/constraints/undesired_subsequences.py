import regex as re


class UndesiredSubSequenceFinder(object):
    def __init__(self, file):
        if file == []:
            file = None
        self.file = file
        if file is not None:
            with open(file) as f:
                list_of_subsequences = [seq.strip() for seq in f.readlines()]
        else:
            list_of_subsequences = []
        self.pattern_obj = re.compile("|".join(list_of_subsequences))

    def undesired_subsequences_val(self, sequence):
        res = [0.0] * len(sequence)
        for m in self.pattern_obj.finditer(sequence, overlapped=True):
            res[m.span()[0]:m.span()[1]] = [1.0] * (m.span()[1] - m.span()[0])
        return res

    def undesired_subsequences(self, sequence):
        return max(self.undesired_subsequences_val(sequence)) > 0.0


if __name__ == "__main__":
    uFinder = UndesiredSubSequenceFinder("undesired_subsequences.txt")
    print(uFinder.undesired_subsequences("ACCCAGTGACCCTATAGCCCCGTAGCTGCCACCCACCCACCCCAAA"))
    print(uFinder.undesired_subsequences_val("ACCCAGTGACCCTATAGCCCCGTAGCTGCCACCCACCCACCCCAAA"))
