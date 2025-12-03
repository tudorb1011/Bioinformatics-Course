import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path


def read_fasta_file(filename):
    """Read sequences from FASTA file."""
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()  # This removes \r\n or \n
            if not line:
                continue
                
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                
                # Start new sequence - use full header as ID
                current_id = line[1:].strip()
                current_seq = []
            else:
                # Append sequence line
                current_seq.append(line.upper())
        
        # Save last sequence
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences


def calculate_cg_content(sequence):
    """Calculate C+G percentage."""
    if len(sequence) == 0:
        return 0.0
    c_count = sequence.count('C')
    g_count = sequence.count('G')
    return ((c_count + g_count) / len(sequence)) * 100


def calculate_kappa_ic(sequence):
    """Calculate Kappa Index of Coincidence."""
    if len(sequence) <= 1:
        return 0.0
    
    N = len(sequence)
    nucleotides = ['A', 'C', 'G', 'T']
    counts = {nuc: sequence.count(nuc) for nuc in nucleotides}
    
    numerator = sum(count * (count - 1) for count in counts.values())
    denominator = N * (N - 1)
    
    return (numerator / denominator) * 85


def sliding_window_analysis(sequence, window_size=30):
    """Perform sliding window analysis."""
    cg_values = []
    ic_values = []
    positions = []
    
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        
        cg = calculate_cg_content(window)
        ic = calculate_kappa_ic(window)
        
        cg_values.append(cg)
        ic_values.append(ic)
        positions.append(i + window_size // 2)
    
    return cg_values, ic_values, positions


def calculate_center_of_weight(cg_values, ic_values):
    """Calculate center of weight (centroid)."""
    if len(cg_values) == 0:
        return 0.0, 0.0
    
    weights = np.array([cg + ic for cg, ic in zip(cg_values, ic_values)])
    
    if weights.sum() == 0:
        return np.mean(cg_values), np.mean(ic_values)
    
    weights = weights / weights.sum()
    center_cg = np.dot(weights, cg_values)
    center_ic = np.dot(weights, ic_values)
    
    return center_cg, center_ic


def generate_digital_stain(seq_id, sequence, window_size, output_dir):
    """Generate and save digital stain for a promoter sequence."""
    
    # Perform analysis
    cg_values, ic_values, positions = sliding_window_analysis(sequence, window_size)
    
    if len(cg_values) == 0:
        return None
    
    # Calculate center
    center_cg, center_ic = calculate_center_of_weight(cg_values, ic_values)
    
    # Create plot
    plt.figure(figsize=(10, 8))
    
    scatter = plt.scatter(cg_values, ic_values, c=positions, cmap='viridis', 
                         s=50, alpha=0.6, edgecolors='black', linewidth=0.5)
    
    plt.colorbar(scatter, label='Position in Sequence')
    
    plt.plot(center_cg, center_ic, 'r*', markersize=20, 
            label=f'Center (CG={center_cg:.2f}, IC={center_ic:.2f})',
            markeredgecolor='black', markeredgewidth=1.5)
    
    plt.xlabel('C+G% Content', fontsize=12)
    plt.ylabel('Kappa Index of Coincidence', fontsize=12)
    plt.title(f'Digital Stain - {seq_id}', fontsize=14, fontweight='bold')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save plot
    output_file = os.path.join(output_dir, f'{seq_id}_digital_stain.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    # Return statistics
    overall_cg = calculate_cg_content(sequence)
    overall_ic = calculate_kappa_ic(sequence)
    
    return {
        'id': seq_id,
        'length': len(sequence),
        'overall_cg': overall_cg,
        'overall_ic': overall_ic,
        'center_cg': center_cg,
        'center_ic': center_ic
    }


def create_summary_plot(results, output_dir):
    """Create a summary plot showing all promoter centers."""
    
    plt.figure(figsize=(12, 10))
    
    # Extract data
    centers_cg = [r['center_cg'] for r in results]
    centers_ic = [r['center_ic'] for r in results]
    overall_cg = [r['overall_cg'] for r in results]
    
    # Create scatter plot colored by overall CG content
    scatter = plt.scatter(centers_cg, centers_ic, c=overall_cg, 
                         cmap='RdYlGn', s=80, alpha=0.6, 
                         edgecolors='black', linewidth=0.5)
    
    cbar = plt.colorbar(scatter)
    cbar.set_label('Overall C+G% Content', rotation=270, labelpad=20)
    
    plt.xlabel('C+G% Content (Center of Weight)', fontsize=12)
    plt.ylabel('Kappa Index of Coincidence (Center of Weight)', fontsize=12)
    plt.title(f'Digital Stain Centers - All Promoters (n={len(results)})', 
             fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    output_file = os.path.join(output_dir, 'all_promoters_summary.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {output_file}")


def main():
    # Configuration
    input_file = "promotori_lista_completa.txt"
    window_size = 30
    output_dir = "digital_stains"
    max_individual_plots = 10  # Only save first 10 individual images
    
    print("="*70)
    print("DNA PROMOTER PATTERN ANALYZER")
    print("="*70)
    
    # Check input file
    if not os.path.exists(input_file):
        print(f"\nError: File '{input_file}' not found!")
        print("Please place the promoter file in the same directory.")
        return
    
    # Create output directory
    Path(output_dir).mkdir(exist_ok=True)
    
    # Read sequences
    print(f"\nReading: {input_file}")
    sequences = read_fasta_file(input_file)
    print(f"Found: {len(sequences)} promoter sequences")
    print(f"Note: Will save only first {max_individual_plots} individual plots")
    print(f"      All sequences will be included in summary plot")
    
    # Process sequences
    print(f"\nProcessing sequences...")
    print("Window size: 30 bp\n")
    
    results = []
    for idx, (seq_id, sequence) in enumerate(sequences.items(), 1):
        # Only save individual plots for first N sequences
        save_individual = (idx <= max_individual_plots)
        
        if save_individual:
            result = generate_digital_stain(seq_id, sequence, window_size, output_dir)
        else:
            # Just compute statistics without saving plot
            cg_values, ic_values, positions = sliding_window_analysis(sequence, window_size)
            if len(cg_values) > 0:
                center_cg, center_ic = calculate_center_of_weight(cg_values, ic_values)
                result = {
                    'id': seq_id,
                    'length': len(sequence),
                    'overall_cg': calculate_cg_content(sequence),
                    'overall_ic': calculate_kappa_ic(sequence),
                    'center_cg': center_cg,
                    'center_ic': center_ic
                }
            else:
                result = None
        
        if result:
            results.append(result)
            
            # Progress updates
            if idx <= max_individual_plots:
                print(f"  [{idx}/{len(sequences)}] {seq_id}: CG={result['overall_cg']:.2f}%, IC={result['overall_ic']:.2f} (plot saved)")
            elif idx % 100 == 0:
                print(f"  [{idx}/{len(sequences)}] Processed...")
    
    # Create summary plot with ALL sequences
    print(f"\nGenerating summary plot with all {len(results)} sequences...")
    create_summary_plot(results, output_dir)
    
    print(f"\n{'='*70}")
    print(f"Complete! Processed {len(results)} promoter sequences")
    print(f"\nOutput folder: {output_dir}/")
    print(f"  - Individual plots: {min(max_individual_plots, len(results))} files")
    print(f"  - Summary plot: all_promoters_summary.png")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
