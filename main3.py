# using ai to adapt ypur current algorithms from main1 and main2, in order to make an app that takes a FASTA file and
# read the seq content from it and display the r procentage for the symbols present in the al[phabet of the seq. note: fasta is
# a file format that contains DNA, ARN or proteins seq. Thus it contains the information for your input.

def read_fasta(file_path: str) -> str:
    """Read a FASTA file and return the sequence content as a single string."""
    sequence = []
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):  # ignoră header-ul și liniile goale
                continue
            sequence.append(line.upper())  # păstrăm literele mari (DNA/ARN/Proteine)
    return "".join(sequence)


def calculate_percentage(sequence: str) -> dict:
    """Calculate frequency and percentage of each symbol in the sequence."""
    alphabet = sorted(set(sequence))
    total_length = len(sequence)
    percentage_dict = {}

    for letter in alphabet:
        count = sequence.count(letter)
        percentage = (count / total_length) * 100
        percentage_dict[letter] = {
            'count': count,
            'percentage': round(percentage, 2)
        }

    return percentage_dict


if __name__ == "__main__":
    fasta_file = "example.fasta"  # aici pui numele fișierului FASTA
    seq = read_fasta(fasta_file)
    result = calculate_percentage(seq)

    print("Alphabet:", sorted(set(seq)))
    print("Total length:", len(seq))
    print("\nComposition:")
    for letter, stats in result.items():
        print(f"Letter: {letter}, Count: {stats['count']}, Percentage: {stats['percentage']}%")
