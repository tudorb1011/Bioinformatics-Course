
def find_alphabet(sequence: str) -> list:

    return sorted(set(sequence))


if __name__ == "__main__":
    seq1 = "abbbbac"
    seq2 = "adcbbbac"
    print(find_alphabet(seq1))
    print(find_alphabet(seq2))