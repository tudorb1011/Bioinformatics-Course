from genetic_translator import translate_sequence, format_output, dna_to_rna

def test_translations():
    """Test the translator with example sequences"""

    # Test cases
    test_cases = [
        {
            'name': 'Simple DNA sequence',
            'sequence': 'ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG',
            'description': 'DNA sequence starting with ATG (Met)'
        },
        {
            'name': 'RNA sequence',
            'sequence': 'AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG',
            'description': 'RNA sequence starting with AUG (Met)'
        },
        {
            'name': 'DNA without start codon',
            'sequence': 'GCCATTGTAATGGGCCGCTGA',
            'description': 'DNA sequence without start codon'
        },
        {
            'name': 'Sequence with stop codon',
            'sequence': 'ATGGCCTAAGTAATG',
            'description': 'Sequence ending with UAA stop codon'
        }
    ]
    
    for i, test in enumerate(test_cases, 1):
        print(f"\n{'=' * 70}")
        print(f"Test {i}: {test['name']}")
        print(f"Description: {test['description']}")
        print(f"{'=' * 70}")
        
        sequence = test['sequence']
        print(f"\nOriginal Sequence: {sequence}")
        
        # Convert to RNA for display
        rna = dna_to_rna(sequence)
        print(f"RNA Sequence:      {rna}")
        
        # Translate with start codon
        amino_acids_with_start = translate_sequence(sequence, start_codon='AUG')
        print(f"\n--- Translation (with start codon search) ---")
        if amino_acids_with_start:
            print(f"Amino acids (3-letter): {format_output(amino_acids_with_start, 'full')}")
            print(f"Amino acids (1-letter): {format_output(amino_acids_with_start, 'single')}")
            print(f"Count: {len(amino_acids_with_start)} amino acids")
        else:
            print("No amino acids produced (start codon not found or sequence too short)")
        
        # Translate without start codon requirement
        amino_acids_no_start = translate_sequence(sequence, start_codon=None)
        print(f"\n--- Translation (from beginning, no start codon) ---")
        if amino_acids_no_start:
            print(f"Amino acids (3-letter): {format_output(amino_acids_no_start, 'full')}")
            print(f"Amino acids (1-letter): {format_output(amino_acids_no_start, 'single')}")
            print(f"Count: {len(amino_acids_no_start)} amino acids")
        else:
            print("No amino acids produced")
    
    print(f"\n{'=' * 70}")
    print("All tests completed!")
    print(f"{'=' * 70}\n")


if __name__ == "__main__":
    test_translations()
