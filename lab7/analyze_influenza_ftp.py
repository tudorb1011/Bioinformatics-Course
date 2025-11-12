import matplotlib.pyplot as plt
from collections import Counter
import os

# Check if influenza.fna exists
if not os.path.exists('influenza.fna'):
    print("ERROR: influenza.fna not found!")
    print("\nPlease download it first:")
    print("  wget https://ftp.ncbi.nlm.nih.gov/genomes/INFLUENZA/influenza.fna")
    print("  or")
    print("  curl -O https://ftp.ncbi.nlm.nih.gov/genomes/INFLUENZA/influenza.fna")
    exit(1)

print("Parsing influenza.fna...")
genomes = []
current_seq = []
current_header = None

with open("influenza.fna", 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if current_header and current_seq:
                seq = ''.join(current_seq).upper()
                # Only keep valid DNA sequences
                seq = ''.join(c for c in seq if c in 'ATGC')
                if 1000 <= len(seq) <= 3000:
                    genomes.append({'header': current_header, 'seq': seq})
                    if len(genomes) == 10:
                        break
            current_header = line
            current_seq = []
        else:
            current_seq.append(line)

print(f"Found {len(genomes)} genomes (1000-3000 bp)")

def analyze_repetitions(sequence):
    repetitions = Counter()
    for length in range(6, 11):
        for i in range(len(sequence) - length + 1):
            pattern = sequence[i:i + length]
            repetitions[pattern] += 1
    return {p: c for p, c in repetitions.items() if c > 1}

def plot_genome(repetitions, title, filename, seq_length):
    top_20 = sorted(repetitions.items(), key=lambda x: x[1], reverse=True)[:20]
    patterns, frequencies = zip(*top_20)
    
    length_freq = Counter()
    for p, c in repetitions.items():
        length_freq[len(p)] += c
    lengths = sorted(length_freq.keys())
    counts = [length_freq[l] for l in lengths]
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
    
    ax1.bar(range(len(patterns)), frequencies, color='steelblue', alpha=0.8)
    ax1.set_xlabel('DNA Pattern', fontsize=11)
    ax1.set_ylabel('Frequency', fontsize=11)
    ax1.set_title(f'Top 20 Repetitions - {title}\n{seq_length} bp', fontsize=11, fontweight='bold')
    ax1.set_xticks(range(len(patterns)))
    ax1.set_xticklabels(patterns, rotation=45, ha='right', fontsize=8)
    ax1.grid(axis='y', alpha=0.3)
    for i, v in enumerate(frequencies):
        ax1.text(i, v, str(v), ha='center', va='bottom', fontsize=8)
    
    ax2.bar(lengths, counts, color='coral', alpha=0.8)
    ax2.set_xlabel('Pattern Length (bp)', fontsize=11)
    ax2.set_ylabel('Total Frequency', fontsize=11)
    ax2.set_title('Frequency by Length', fontsize=11, fontweight='bold')
    ax2.set_xticks(lengths)
    ax2.grid(axis='y', alpha=0.3)
    for l, c in zip(lengths, counts):
        ax2.text(l, c, str(c), ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

def create_comparison(all_results):
    fig, axes = plt.subplots(2, 5, figsize=(20, 8))
    fig.suptitle('10 Influenza Genomes - Repetition Frequency Comparison', fontsize=16, fontweight='bold')
    
    for idx, data in enumerate(all_results):
        ax = axes[idx // 5, idx % 5]
        repetitions = data['repetitions']
        
        length_freq = Counter()
        for p, c in repetitions.items():
            length_freq[len(p)] += c
        lengths = sorted(length_freq.keys())
        counts = [length_freq[l] for l in lengths]
        
        ax.bar(lengths, counts, color='steelblue', alpha=0.7)
        ax.set_title(f"Genome {idx+1}\n{data['length']} bp", fontsize=9)
        ax.set_xlabel('Length (bp)', fontsize=8)
        ax.set_ylabel('Frequency', fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('influenza_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()

os.makedirs('influenza_genomes', exist_ok=True)

print("\n" + "="*70)
print("ANALYZING 10 INFLUENZA GENOMES FROM NCBI FTP")
print("="*70)

all_results = []

for i, genome in enumerate(genomes, 1):
    print(f"\n[{i}/10] Processing genome {i}")
    
    header_short = genome['header'][:60]
    
    # Save FASTA
    with open(f"influenza_genomes/genome_{i}.fasta", 'w') as f:
        f.write(f"{genome['header']}\n{genome['seq']}\n")
    
    # Analyze
    repetitions = analyze_repetitions(genome['seq'])
    
    all_results.append({
        'header': header_short,
        'length': len(genome['seq']),
        'repetitions': repetitions
    })
    
    # Plot
    plot_file = f"influenza_genomes/genome_{i}_repetitions.png"
    plot_genome(repetitions, f"Genome {i}", plot_file, len(genome['seq']))
    
    print(f"  Length: {len(genome['seq'])} bp")
    print(f"  Patterns: {len(repetitions)}")
    print(f"  Saved: {plot_file}")

print("\n" + "="*70)
print("Creating comparison plot...")
create_comparison(all_results)
print("Saved: influenza_comparison.png")

print("\n" + "="*70)
print("SUMMARY")
print("="*70)
for i, data in enumerate(all_results, 1):
    print(f"Genome {i}: {data['length']} bp, {len(data['repetitions'])} unique patterns")
    print(f"  {data['header']}")

print("\n" + "="*70)
print("COMPLETE!")
print("="*70)
