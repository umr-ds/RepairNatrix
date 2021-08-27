def homopolymer(sequence, count):
    prev_char = sequence[0]
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


if __name__ == "__main__":
    print(homopolymer("AAAAAAAAAACACACTTTTTTTTAAAAAAAAAAA",5))
    print(homopolymer("AAAAGGGGGGCGGGGAGCCCTTTTTCGCGCCCCCCGGGGTTTTTT",5))
