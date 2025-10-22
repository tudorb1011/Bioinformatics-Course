import re
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Genetic code dictionary
GENETIC_CODE = {
    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop', 'UGG': 'Trp',
    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
}


def read_fasta(filename):
    """Read FASTA file and return sequence"""
    sequence = ""
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('>'):
                sequence += line.upper()
    return sequence


def dna_to_rna(sequence):
    """Convert DNA to RNA"""
    return sequence.replace('T', 'U')


def extract_codons(sequence):
    """Extract all codons from sequence"""
    # Convert to RNA if needed
    if 'T' in sequence:
        sequence = dna_to_rna(sequence)
    
    # Extract codons (groups of 3)
    codons = []
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if len(codon) == 3 and codon in GENETIC_CODE:
            codons.append(codon)
    
    return codons


def count_codons(codons):
    """Count codon frequencies"""
    return Counter(codons)


def count_amino_acids(codons):
    """Count amino acid frequencies from codons"""
    amino_acids = []
    for codon in codons:
        if codon in GENETIC_CODE:
            aa = GENETIC_CODE[codon]
            if aa != 'Stop':
                amino_acids.append(aa)
    return Counter(amino_acids)


def plot_top_codons(codon_counts, title, filename, top_n=10):
    """Plot top N most frequent codons"""
    # Get top N codons
    top_codons = codon_counts.most_common(top_n)
    
    codons = [codon for codon, count in top_codons]
    counts = [count for codon, count in top_codons]
    
    # Add amino acid labels
    labels = [f"{codon}\n({GENETIC_CODE.get(codon, '?')})" for codon in codons]
    
    # Create plot
    plt.figure(figsize=(12, 6))
    bars = plt.bar(range(len(codons)), counts, color='steelblue', edgecolor='black')
    
    # Add value labels on bars
    for i, (bar, count) in enumerate(zip(bars, counts)):
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                f'{count:,}',
                ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    plt.xlabel('Codon (Amino Acid)', fontsize=12, fontweight='bold')
    plt.ylabel('Frequency', fontsize=12, fontweight='bold')
    plt.title(title, fontsize=14, fontweight='bold')
    plt.xticks(range(len(codons)), labels, fontsize=10)
    plt.grid(axis='y', alpha=0.3, linestyle='--')
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"✓ Chart saved: {filename}")
    plt.close()


def plot_comparison(covid_counts, flu_counts, filename, top_n=15):
    """Plot comparison of top codons between COVID-19 and Influenza"""
    # Get all codons from both
    all_codons = set(list(dict(covid_counts.most_common(top_n)).keys()) + 
                     list(dict(flu_counts.most_common(top_n)).keys()))
    
    # Get counts for all codons
    codons_list = sorted(all_codons)
    covid_values = [covid_counts.get(codon, 0) for codon in codons_list]
    flu_values = [flu_counts.get(codon, 0) for codon in codons_list]
    
    # Normalize to percentages
    total_covid = sum(covid_counts.values())
    total_flu = sum(flu_counts.values())
    covid_percentages = [(count / total_covid) * 100 for count in covid_values]
    flu_percentages = [(count / total_flu) * 100 for count in flu_values]
    
    # Create plot
    fig, ax = plt.subplots(figsize=(14, 8))
    
    x = range(len(codons_list))
    width = 0.35
    
    bars1 = ax.bar([i - width/2 for i in x], covid_percentages, width, 
                    label='COVID-19', color='#e74c3c', alpha=0.8, edgecolor='black')
    bars2 = ax.bar([i + width/2 for i in x], flu_percentages, width, 
                    label='Influenza', color='#3498db', alpha=0.8, edgecolor='black')
    
    # Add labels
    labels = [f"{codon}\n({GENETIC_CODE.get(codon, '?')})" for codon in codons_list]
    ax.set_xlabel('Codon (Amino Acid)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Frequency (%)', fontsize=12, fontweight='bold')
    ax.set_title('COVID-19 vs Influenza: Most Frequent Codons Comparison', 
                 fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.legend(fontsize=11)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"✓ Comparison chart saved: {filename}")
    plt.close()


def print_statistics(name, sequence, codon_counts, aa_counts):
    """Print statistics for a genome"""
    print(f"\n{'='*70}")
    print(f"{name} STATISTICS")
    print(f"{'='*70}")
    print(f"Genome length: {len(sequence):,} nucleotides")
    print(f"Total codons: {sum(codon_counts.values()):,}")
    print(f"Unique codons: {len(codon_counts)}")
    print(f"Total amino acids: {sum(aa_counts.values()):,}")
    
    # Top 10 codons
    print(f"\n--- Top 10 Most Frequent Codons ---")
    for i, (codon, count) in enumerate(codon_counts.most_common(10), 1):
        aa = GENETIC_CODE.get(codon, '?')
        percentage = (count / sum(codon_counts.values())) * 100
        print(f"{i:2d}. {codon} ({aa:3s}): {count:5,} ({percentage:5.2f}%)")
    
    # Top 3 amino acids
    print(f"\n--- Top 3 Most Frequent Amino Acids ---")
    for i, (aa, count) in enumerate(aa_counts.most_common(3), 1):
        percentage = (count / sum(aa_counts.values())) * 100
        print(f"{i}. {aa}: {count:,} ({percentage:.2f}%)")


def compare_codons(covid_counts, flu_counts):
    """Compare and find common frequent codons"""
    print(f"\n{'='*70}")
    print("COMPARISON: MOST FREQUENT CODONS IN BOTH GENOMES")
    print(f"{'='*70}")
    
    # Get top 20 from each
    covid_top = set([codon for codon, count in covid_counts.most_common(20)])
    flu_top = set([codon for codon, count in flu_counts.most_common(20)])
    
    # Find common codons
    common_codons = covid_top & flu_top
    
    print(f"\nCommon codons in top 20 of both genomes: {len(common_codons)}")
    
    # Sort by combined frequency
    combined_freq = []
    for codon in common_codons:
        covid_freq = covid_counts.get(codon, 0)
        flu_freq = flu_counts.get(codon, 0)
        combined_freq.append((codon, covid_freq, flu_freq, covid_freq + flu_freq))
    
    combined_freq.sort(key=lambda x: x[3], reverse=True)
    
    print(f"\n{'Codon':6s} {'AA':4s} {'COVID-19':>12s} {'Influenza':>12s} {'Combined':>12s}")
    print('-' * 60)
    
    for codon, covid_freq, flu_freq, combined in combined_freq[:15]:
        aa = GENETIC_CODE.get(codon, '?')
        print(f"{codon:6s} {aa:4s} {covid_freq:12,} {flu_freq:12,} {combined:12,}")


def main():
    """Main analysis function"""
    print("="*70)
    print("CODON FREQUENCY ANALYSIS: COVID-19 vs INFLUENZA")
    print("="*70)
    
    # File names (expected in working directory)
    covid_file = "covid.fasta"
    flu_file = "influenzaD.fasta"
    
    try:
        # Read sequences
        print("\n[1/7] Reading COVID-19 genome...")
        covid_seq = read_fasta(covid_file)
        print(f"✓ Loaded: {len(covid_seq):,} nucleotides")
        
        print("\n[2/7] Reading Influenza genome...")
        flu_seq = read_fasta(flu_file)
        print(f"✓ Loaded: {len(flu_seq):,} nucleotides")
        
        # Extract codons
        print("\n[3/7] Extracting codons...")
        covid_codons = extract_codons(covid_seq)
        flu_codons = extract_codons(flu_seq)
        print(f"✓ COVID-19 codons: {len(covid_codons):,}")
        print(f"✓ Influenza codons: {len(flu_codons):,}")
        
        # Count frequencies
        print("\n[4/7] Counting codon frequencies...")
        covid_codon_counts = count_codons(covid_codons)
        flu_codon_counts = count_codons(flu_codons)
        print("✓ Frequencies calculated")
        
        # Count amino acids
        print("\n[5/7] Counting amino acid frequencies...")
        covid_aa_counts = count_amino_acids(covid_codons)
        flu_aa_counts = count_amino_acids(flu_codons)
        print("✓ Amino acids counted")
        
        # Print statistics
        print_statistics("COVID-19", covid_seq, covid_codon_counts, covid_aa_counts)
        print_statistics("INFLUENZA", flu_seq, flu_codon_counts, flu_aa_counts)
        
        # Create charts
        print(f"\n[6/7] Creating charts...")
        
        # A) COVID-19 top 10
        plot_top_codons(covid_codon_counts, 
                       "COVID-19: Top 10 Most Frequent Codons",
                       "covid19_top10_codons.png", top_n=10)
        
        # B) Influenza top 10
        plot_top_codons(flu_codon_counts,
                       "Influenza: Top 10 Most Frequent Codons",
                       "influenza_top10_codons.png", top_n=10)
        
        # C) Comparison
        plot_comparison(covid_codon_counts, flu_codon_counts,
                       "covid_vs_influenza_comparison.png", top_n=15)
        
        # D) Compare results
        print(f"\n[7/7] Comparing genomes...")
        compare_codons(covid_codon_counts, flu_codon_counts)
        
        print(f"\n{'='*70}")
        print("ANALYSIS COMPLETE!")
        print(f"{'='*70}")
        print("\nGenerated files:")
        print("  - covid19_top10_codons.png")
        print("  - influenza_top10_codons.png")
        print("  - covid_vs_influenza_comparison.png")
        print()
        
    except FileNotFoundError as e:
        print(f"\n❌ ERROR: File not found - {e}")
        print("\nPlease make sure these files are in the working directory:")
        print("  - covid.fasta")
        print("  - influenzaD.fasta")
        print()
    except Exception as e:
        print(f"\n❌ ERROR: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
