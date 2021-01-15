from scripts import error_sources


def annotate(sequence):
    tmp = dict()
    for method in error_sources.active_rules:
        tmp[method.__name__] = method(sequence)
    return tmp


if __name__ == "__main__":
    print(annotate("ACCGCACGTGTGGCGCCGCCCCACTAGCTAGCACCCCCCCCCCCCCCACTGACCCCTAGCTTGTATAATGAGCTAGCTGCTACGCTCGCGCCGGCC"))
