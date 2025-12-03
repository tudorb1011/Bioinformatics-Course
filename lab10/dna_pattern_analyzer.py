import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple


def calculate_cg_content(sequence: str) -> float:

    if len(sequence) == 0:
        return 0.0

    c_count = sequence.upper().count('C')
    g_count = sequence.upper().count('G')
    cg_percentage = ((c_count + g_count) / len(sequence)) * 100

    return cg_percentage


def calculate_index_of_coincidence(sequence: str) -> float:

    if len(sequence) <= 1:
        return 0.0

    sequence = sequence.upper()
    N = len(sequence)

    # Count occurrences of each nucleotide
    nucleotides = ['A', 'C', 'G', 'T']
    counts = {nuc: sequence.count(nuc) for nuc in nucleotides}

    # Calculate IC using standard formula: sum(ni * (ni-1)) / (N * (N-1))
    numerator = sum(count * (count - 1) for count in counts.values())
    denominator = N * (N - 1)

    # Scale by 85 to match PromKappa output
    kappa = (numerator / denominator) * 85

    return kappa


def sliding_window_analysis(sequence: str, window_size: int) -> Tuple[List[float], List[float], List[int]]:

    cg_values = []
    ic_values = []
    positions = []

    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]

        cg = calculate_cg_content(window)
        ic = calculate_index_of_coincidence(window)

        cg_values.append(cg)
        ic_values.append(ic)
        positions.append(i + window_size // 2)  # Center position of window

    return cg_values, ic_values, positions


def calculate_center_of_weight(cg_values: List[float], ic_values: List[float],
                               positions: List[int]) -> Tuple[float, float]:

    if len(cg_values) == 0:
        return 0.0, 0.0

    # Weight by the magnitude of the values
    weights = np.array([cg + ic for cg, ic in zip(cg_values, ic_values)])

    if weights.sum() == 0:
        center_cg = np.mean(cg_values)
        center_ic = np.mean(ic_values)
    else:
        weights = weights / weights.sum()
        center_cg = np.dot(weights, cg_values)
        center_ic = np.dot(weights, ic_values)

    return center_cg, center_ic


def plot_dna_pattern(cg_values: List[float], ic_values: List[float],
                     positions: List[int], title: str = "DNA Pattern Analysis",
                     center: Tuple[float, float] = None):

    plt.figure(figsize=(10, 8))

    # Create scatter plot with color gradient based on position
    scatter = plt.scatter(cg_values, ic_values, c=positions, cmap='viridis',
                          s=50, alpha=0.6, edgecolors='black', linewidth=0.5)

    # Add colorbar to show position along sequence
    cbar = plt.colorbar(scatter)
    cbar.set_label('Position in Sequence', rotation=270, labelpad=20)

    # Plot center of weight if provided
    if center is not None:
        plt.plot(center[0], center[1], 'r*', markersize=20,
                 label=f'Center of Weight\n(CG={center[0]:.2f}, IC={center[1]:.2f})',
                 markeredgecolor='black', markeredgewidth=1.5)
        plt.legend()

    plt.xlabel('C+G% Content', fontsize=12)
    plt.ylabel('Kappa Index of Coincidence', fontsize=12)
    plt.title(title, fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()


def plot_pattern_centers(centers: List[Tuple[float, float, str]]):

    plt.figure(figsize=(10, 8))

    for i, (cg_center, ic_center, label) in enumerate(centers):
        plt.plot(cg_center, ic_center, 'o', markersize=12,
                 label=label, markeredgecolor='black', markeredgewidth=1.5)

    plt.xlabel('C+G% Content (Center)', fontsize=12)
    plt.ylabel('Kappa Index of Coincidence (Center)', fontsize=12)
    plt.title('Pattern Centers Comparison', fontsize=14, fontweight='bold')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()


def main():

    print("=" * 70)
    print("DNA PATTERN ANALYZER")
    print("=" * 70)

    # Step 1: Define the test sequence
    sequence = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"
    print(f"\nSequence length: {len(sequence)} bp")
    print(f"Sequence: {sequence[:50]}..." if len(sequence) > 50 else f"Sequence: {sequence}")

    # Step 2: Define window size
    window_size = 30
    print(f"\nWindow size: {window_size} bp")

    # Step 3 & 4: Perform sliding window analysis
    print("\nPerforming sliding window analysis...")
    cg_values, ic_values, positions = sliding_window_analysis(sequence, window_size)

    # Validation checks
    print(f"\nNumber of windows analyzed: {len(cg_values)}")

    # Calculate overall sequence statistics for validation
    overall_cg = calculate_cg_content(sequence)
    overall_ic = calculate_index_of_coincidence(sequence)

    print(f"\nValidation Results:")
    print(f"  Overall Sequence C+G%: {overall_cg:.2f} (Expected: 29.27)")
    print(f"  Overall Sequence IC:   {overall_ic:.2f} (Expected: 27.53)")
    print(f"\n  Average Window C+G%: {np.mean(cg_values):.2f}")
    print(f"  Average Window IC:   {np.mean(ic_values):.2f}")

    # Show range of values
    print(f"\nValue Ranges:")
    print(f"  C+G%: {min(cg_values):.2f} - {max(cg_values):.2f}")
    print(f"  IC:   {min(ic_values):.2f} - {max(ic_values):.2f}")

    # Step 6: Calculate center of weight
    center_cg, center_ic = calculate_center_of_weight(cg_values, ic_values, positions)
    print(f"\nCenter of Weight:")
    print(f"  C+G%: {center_cg:.2f}")
    print(f"  IC:   {center_ic:.2f}")

    # Step 5: Plot the pattern
    print("\nGenerating pattern plot...")
    plot_dna_pattern(cg_values, ic_values, positions,
                     title="DNA Pattern - Test Sequence",
                     center=(center_cg, center_ic))
    plt.savefig('dna_pattern_main.png', dpi=300, bbox_inches='tight')
    print("  Saved: dna_pattern_main.png")

    # Step 7: Create a chart with pattern centers
    centers = [(center_cg, center_ic, "Test Sequence")]

    # You can add more sequences here for comparison
    # Example: analyze another sequence
    print("\n" + "=" * 70)
    print("Additional Analysis Example")
    print("=" * 70)

    # Example with a different sequence (high CG content)
    high_cg_seq = "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG"
    print(f"\nHigh CG Sequence: {high_cg_seq[:50]}...")
    cg_vals_hcg, ic_vals_hcg, pos_hcg = sliding_window_analysis(high_cg_seq, window_size)
    center_cg_hcg, center_ic_hcg = calculate_center_of_weight(cg_vals_hcg, ic_vals_hcg, pos_hcg)
    print(f"Center: C+G%={center_cg_hcg:.2f}, IC={center_ic_hcg:.2f}")
    centers.append((center_cg_hcg, center_ic_hcg, "High CG Content"))

    # Example with AT-rich sequence
    at_rich_seq = "ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT"
    print(f"\nAT-Rich Sequence: {at_rich_seq[:50]}...")
    cg_vals_at, ic_vals_at, pos_at = sliding_window_analysis(at_rich_seq, window_size)
    center_cg_at, center_ic_at = calculate_center_of_weight(cg_vals_at, ic_vals_at, pos_at)
    print(f"Center: C+G%={center_cg_at:.2f}, IC={center_ic_at:.2f}")
    centers.append((center_cg_at, center_ic_at, "AT-Rich"))

    # Step 7: Plot pattern centers
    print("\nGenerating pattern centers comparison plot...")
    plot_pattern_centers(centers)
    plt.savefig('dna_pattern_centers.png', dpi=300, bbox_inches='tight')
    print("  Saved: dna_pattern_centers.png")

    plt.show()

    print("\n" + "=" * 70)
    print("Analysis Complete!")
    print("=" * 70)
    print("\nTo analyze your own promoter sequence:")
    print("  1. Replace the 'sequence' variable with your promoter sequence")
    print("  2. Adjust window_size if needed")
    print("  3. Run the script again")
    print("\nOutput files saved in current directory")


if __name__ == "__main__":
    main()
