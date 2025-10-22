"""
Genetic Code Translator
Converts DNA/RNA coding sequences into amino acid sequences
"""

# Genetic code dictionary based on the provided table
GENETIC_CODE = {
    # U (First letter)
    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop', 'UGG': 'Trp',
    
    # C (First letter)
    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    
    # A (First letter)
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    
    # G (First letter)
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
}


def dna_to_rna(dna_sequence):
    """Convert DNA sequence to RNA by replacing T with U"""
    return dna_sequence.upper().replace('T', 'U')


def translate_sequence(sequence, start_codon='AUG'):
    """
    Translate a genetic sequence into amino acids
    
    Args:
        sequence: DNA or RNA sequence string
        start_codon: Starting codon for translation (default: AUG for Met)
    
    Returns:
        List of amino acids
    """
    # Convert to RNA if it's DNA
    if 'T' in sequence.upper():
        sequence = dna_to_rna(sequence)
    else:
        sequence = sequence.upper()
    
    # Find start codon if specified
    if start_codon:
        start_codon = start_codon.upper().replace('T', 'U')
        if start_codon in sequence:
            start_index = sequence.find(start_codon)
            sequence = sequence[start_index:]
        else:
            print(f"Warning: Start codon {start_codon} not found. Translating from beginning.")
    
    # Translate sequence
    amino_acids = []
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if len(codon) == 3:
            amino_acid = GENETIC_CODE.get(codon, 'Unknown')
            if amino_acid == 'Stop':
                break
            amino_acids.append(amino_acid)
    
    return amino_acids


def format_output(amino_acids, format_type='full'):
    """
    Format amino acid sequence for output
    
    Args:
        amino_acids: List of amino acids
        format_type: 'full' for full names, 'single' for single letter codes
    
    Returns:
        Formatted string
    """
    if format_type == 'full':
        return '-'.join(amino_acids)
    elif format_type == 'single':
        # Single letter amino acid codes
        single_letter = {
            'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
            'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
            'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
            'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
        }
        return ''.join([single_letter.get(aa, 'X') for aa in amino_acids])
    else:
        return ' '.join(amino_acids)


def main():
    """Main function to run the translator"""
    print("=" * 60)
    print("Genetic Code Translator")
    print("=" * 60)
    print()
    
    # Get user input
    sequence = input("Enter DNA or RNA sequence: ").strip()
    
    if not sequence:
        print("Error: No sequence provided!")
        return
    
    # Ask if user wants to find start codon
    find_start = input("Find start codon (AUG/ATG)? (y/n, default: y): ").strip().lower()
    start_codon = 'AUG' if find_start != 'n' else None
    
    # Translate sequence
    amino_acids = translate_sequence(sequence, start_codon)
    
    if not amino_acids:
        print("\nNo amino acids produced. Check your sequence.")
        return
    
    # Display results
    print("\n" + "=" * 60)
    print("RESULTS")
    print("=" * 60)
    
    # Convert input to RNA for display
    rna_seq = sequence.upper().replace('T', 'U')
    if start_codon and start_codon in rna_seq:
        start_idx = rna_seq.find(start_codon)
        rna_seq = rna_seq[start_idx:]
    
    print(f"\nRNA Sequence: {rna_seq[:min(len(rna_seq), 60)]}{'...' if len(rna_seq) > 60 else ''}")
    print(f"Length: {len(rna_seq)} nucleotides")
    print(f"\nNumber of amino acids: {len(amino_acids)}")
    print(f"\nAmino Acid Sequence (3-letter):")
    print(format_output(amino_acids, 'full'))
    print(f"\nAmino Acid Sequence (1-letter):")
    print(format_output(amino_acids, 'single'))
    print()


if __name__ == "__main__":
    main()
