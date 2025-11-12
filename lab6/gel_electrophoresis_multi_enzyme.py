import random
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# Common restriction enzymes with their recognition sites
RESTRICTION_ENZYMES = {
    'EcoRI': 'GAATTC',
    'BamHI': 'GGATCC',
    'HindIII': 'AAGCTT',
    'PstI': 'CTGCAG',
    'SalI': 'GTCGAC'
}


def fetch_dna_sequence(min_length=1000, max_length=3000):
    """Generate a simulated DNA sequence within the specified length range"""
    length = random.randint(min_length, max_length)
    sequence = generate_random_sequence(length)
    print(f"Generated sequence of length: {length} bp")

    # Display the DNA sequence
    print("\nDNA SEQUENCE")
    print(sequence)
    print()

    return sequence


def generate_random_sequence(length):
    """Generate a random DNA sequence with restriction sites"""
    bases = ['A', 'T', 'G', 'C']
    sequence = ''.join(random.choice(bases) for _ in range(length))

    # Insert random restriction sites to ensure enzyme cutting
    num_sites = random.randint(3, 8)
    for _ in range(num_sites):
        enzyme = random.choice(list(RESTRICTION_ENZYMES.keys()))
        site = RESTRICTION_ENZYMES[enzyme]
        insert_pos = random.randint(0, len(sequence) - len(site))
        sequence = sequence[:insert_pos] + site + sequence[insert_pos + len(site):]

    return sequence


def find_restriction_sites(sequence, recognition_site):
    """Find all positions where a restriction enzyme cuts"""
    positions = []
    site_length = len(recognition_site)

    for i in range(len(sequence) - site_length + 1):
        if sequence[i:i + site_length] == recognition_site:
            positions.append(i)

    return positions


def digest_dna(sequence, enzyme_name):
    """
    Simulate DNA digestion with a restriction enzyme
    Returns list of fragment lengths
    """
    recognition_site = RESTRICTION_ENZYMES[enzyme_name]
    cut_positions = find_restriction_sites(sequence, recognition_site)

    if not cut_positions:
        # No cut sites found, return full length
        return [len(sequence)]

    # Add start and end positions
    cut_positions = [0] + cut_positions + [len(sequence)]

    # Calculate fragment lengths
    fragments = []
    for i in range(len(cut_positions) - 1):
        fragment_length = cut_positions[i + 1] - cut_positions[i]
        if fragment_length > 0:
            fragments.append(fragment_length)

    return fragments


def calculate_migration_distance(fragment_length, max_length=3000, gel_length=10):
    """
    Calculate migration distance based on fragment length
    Shorter fragments migrate further (inverse relationship)
    """
    # Logarithmic relationship: shorter fragments travel further
    # Normalize to gel length
    migration = gel_length * (1 - np.log(fragment_length) / np.log(max_length))
    return migration


def visualize_gel_electrophoresis_multi(enzyme_fragments, title="Gel Electrophoresis - Multiple Enzymes"):
    """
    Create a visual representation of gel electrophoresis with multiple enzyme digestions
    Each enzyme gets its own lane
    """
    num_enzymes = len(enzyme_fragments)
    fig, ax = plt.subplots(figsize=(3 + num_enzymes * 1.5, 10), facecolor='black')
    ax.set_facecolor('black')

    # Gel dimensions
    gel_width = 1.2 * num_enzymes
    gel_height = 10

    # Draw the gel box (dark gray/black background)
    gel_box = patches.Rectangle((0.5, 0), gel_width, gel_height,
                                linewidth=0,
                                facecolor='#2a2a2a')
    ax.add_patch(gel_box)

    # Calculate max fragment length across all enzymes for consistent scaling
    all_fragments = []
    for fragments in enzyme_fragments.values():
        all_fragments.extend(fragments)
    max_fragment_length = max(all_fragments) if all_fragments else 3000

    # Add reference markers on the left side
    markers = [3000, 1500, 500]
    for marker_size in markers:
        marker_distance = calculate_migration_distance(marker_size, max_fragment_length, gel_length=gel_height - 1)
        marker_y = gel_height - 0.5 - marker_distance
        ax.text(0.2, marker_y, f'{marker_size} bp',
                ha='right', va='center', fontsize=11,
                fontweight='bold', color='white')
        # Draw reference line
        ax.plot([0.5, 0.5 + gel_width], [marker_y, marker_y],
                'r--', alpha=0.15, linewidth=0.5)

    # Lane settings
    lane_width = gel_width / num_enzymes

    # Draw DNA bands for each enzyme lane
    for idx, (enzyme_name, fragments) in enumerate(enzyme_fragments.items()):
        lane_center = 0.5 + (idx + 0.5) * lane_width
        lane_x_start = lane_center - lane_width * 0.35
        lane_x_width = lane_width * 0.7

        # Draw enzyme label at the top
        ax.text(lane_center, gel_height + 0.3, enzyme_name,
                ha='center', va='bottom', fontsize=10,
                fontweight='bold', color='white')

        # Draw fragment count at bottom
        ax.text(lane_center, -0.3, f'{len(fragments)} frags',
                ha='center', va='top', fontsize=8,
                color='white', style='italic')

        # Draw each fragment band in this lane
        for fragment_length in fragments:
            # Calculate migration distance
            migration_distance = calculate_migration_distance(fragment_length, max_fragment_length,
                                                              gel_length=gel_height - 1)

            # Band position (from top of gel)
            band_y = gel_height - 0.5 - migration_distance

            # Draw the DNA band as a horizontal white bar
            band_height = 0.10
            band = patches.Rectangle((lane_x_start, band_y - band_height / 2),
                                     lane_x_width, band_height,
                                     linewidth=0,
                                     facecolor='white',
                                     alpha=0.9)
            ax.add_patch(band)

            # Add label with bp value on the right side of the band
            ax.text(lane_center + lane_width * 0.45, band_y, f'{fragment_length}',
                    ha='left', va='center', fontsize=8,
                    color='white', fontweight='normal')

    # Formatting
    ax.set_xlim(-0.2, gel_width + 0.8)
    ax.set_ylim(-1, gel_height + 1)
    ax.set_aspect('equal')
    ax.axis('off')

    plt.tight_layout(pad=0.1)
    plt.savefig('gel_electrophoresis_multi.png', dpi=300, bbox_inches='tight', facecolor='black')
    print("\nMulti-enzyme gel electrophoresis visualization saved as 'gel_electrophoresis_multi.png'")

    return fig


def print_enzyme_digest_info(enzyme_fragments, sequence_length):
    """Print information about enzyme digestion results"""
    print("\nRESTRICTION ENZYME DIGESTION ANALYSIS")
    print(f"Original DNA sequence length: {sequence_length} bp\n")

    for enzyme_name, fragments in enzyme_fragments.items():
        recognition_site = RESTRICTION_ENZYMES[enzyme_name]
        num_cuts = len(fragments) - 1

        print(f"{enzyme_name} (Recognition site: {recognition_site})")
        print(f"  Cut sites: {num_cuts}")
        print(f"  Fragments: {len(fragments)}")
        print(f"  Sizes (bp): {sorted(fragments, reverse=True)}")
        print()

    print("SUMMARY - Fragments sorted by size:")
    for enzyme_name, fragments in sorted(enzyme_fragments.items(),
                                         key=lambda x: len(x[1]), reverse=True):
        sorted_frags = sorted(fragments, reverse=True)
        print(f"{enzyme_name:10} ({len(fragments)} fragments): {sorted_frags}")


def main():
    """Main execution function"""
    print("GEL ELECTROPHORESIS SIMULATION - MULTIPLE RESTRICTION ENZYMES")
    print()

    # Generate DNA sequence
    dna_sequence = fetch_dna_sequence(min_length=1000, max_length=3000)

    print(f"Testing with {len(RESTRICTION_ENZYMES)} restriction enzymes:")
    for enzyme, site in RESTRICTION_ENZYMES.items():
        print(f"  {enzyme}: {site}")
    print()

    # Digest DNA with each enzyme
    enzyme_fragments = {}
    for enzyme_name in RESTRICTION_ENZYMES.keys():
        fragments = digest_dna(dna_sequence, enzyme_name)
        enzyme_fragments[enzyme_name] = fragments
        print(f"{enzyme_name}: generated {len(fragments)} fragments")

    # Print detailed analysis
    print_enzyme_digest_info(enzyme_fragments, len(dna_sequence))

    # Visualize gel electrophoresis with all enzymes
    print("\nCreating multi-enzyme gel electrophoresis visualization...")
    visualize_gel_electrophoresis_multi(enzyme_fragments)


if __name__ == "__main__":
    main()