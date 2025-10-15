# a DNA sequence is given: S = 'ACGGGCATATGCGC '. make an app witch is able to show the precentage of the components from the
# alphabet of the sequence S. In other words, the input of the seq s and the output is the alphabet of the seq. and the
# precentage of each letter in the alphabet found in the seq s (freqency and relative freqency)

def calculate_percentage(sequence: str) -> dict:
    alphabet = sorted(set(sequence))
    total_length = len(sequence)
    percentage_dict = {}

    for letter in alphabet:
        count = sequence.count(letter)
        percentage = (count / total_length) * 100  #ca sa fie in procente nu in fractiune
        percentage_dict[letter] = {
            'count': count,
            'percentage': percentage
        }

    return percentage_dict

if __name__ == "__main__":
    seq = "ACGGGCATATGCGC"
    result = calculate_percentage(seq)
    for letter, stats in result.items():
        print("Letter: " + letter + ", Count: " + str(stats['count']) + ", Percentage: " + str(stats['percentage']) + "%")
    print("Alphabet: " + str(sorted(set(seq))))