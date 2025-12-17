import math
import glob
import os
import matplotlib.pyplot as plt

# --- 1. CONFIGURATION & MOTIF DATA ---
MOTIF_LENGTH = 9
BASES = ['A', 'C', 'G', 'T']
BACKGROUND_FREQ = 0.25
GENOME_FOLDER = "genomes"  # <--- Target folder name

motifs = [
    "GAGGTAAAC", "TCCGTAAGT", "CAGGTTGGA", "ACAGTCAGT", "TAGGTCATT",
    "TAGGTACTG", "ATGGTAACT", "CAGGTATAC", "TGTGTGAGT", "AAGGTAAGT"
]


# --- 2. FUNCTIONS ---

def build_log_likelihood_matrix(motifs):
    """Builds the PWM (Log-Likelihood) from the training motifs."""
    num_sequences = len(motifs)
    count_matrix = {b: [0] * MOTIF_LENGTH for b in BASES}
    prob_matrix = {b: [0.0] * MOTIF_LENGTH for b in BASES}
    ll_matrix = {b: [0.0] * MOTIF_LENGTH for b in BASES}

    for seq in motifs:
        for i, char in enumerate(seq):
            if char in count_matrix:
                count_matrix[char][i] += 1

    for b in BASES:
        for i in range(MOTIF_LENGTH):
            prob = count_matrix[b][i] / num_sequences
            prob_matrix[b][i] = prob
            if prob > 0:
                ll_matrix[b][i] = math.log(prob / BACKGROUND_FREQ)
            else:
                ll_matrix[b][i] = -99.0
    return ll_matrix


def read_fasta_file(filepath):
    """Reads a FASTA file from a specific path."""
    sequence = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                # Skip headers and empty lines
                if not line or line.startswith(">"):
                    continue
                sequence.append(line)
        return "".join(sequence).upper()
    except FileNotFoundError:
        print(f"Error: File {filepath} not found.")
        return ""


def scan_genome(sequence, ll_matrix):
    """Scans the sequence and calculates scores."""
    scores = []
    for i in range(len(sequence) - MOTIF_LENGTH + 1):
        window = sequence[i: i + MOTIF_LENGTH]
        current_score = 0
        valid_window = True

        for pos, char in enumerate(window):
            if char not in ll_matrix:
                valid_window = False
                break
            current_score += ll_matrix[char][pos]

        if valid_window:
            scores.append(current_score)
        else:
            scores.append(-99.0)

    return scores


def plot_genome_signal(filepath, scores):
    """Generates and saves a chart for the genome signal."""
    # Extract just the filename for the title (remove 'genomes/' part)
    filename = os.path.basename(filepath)

    plt.figure(figsize=(10, 5))
    plt.plot(scores, color='blue', linewidth=0.5, label='Motif Score')
    plt.axhline(y=0, color='red', linestyle='--', linewidth=1, label='Threshold (0)')

    plt.title(f"Motif Scan Signal: {filename}")
    plt.xlabel("Position in Genome")
    plt.ylabel("Log-Likelihood Score")
    plt.legend()
    plt.grid(True, alpha=0.3)

    # Save chart in the same folder as the genome, or main folder?
    # Let's save it in the main folder to keep 'genomes' clean,
    # or you can change path below.
    output_image = f"scan_result_{filename}.png"
    plt.savefig(output_image)
    print(f"   [Chart saved as {output_image}]")
    plt.close()


# --- 3. MAIN EXECUTION ---

def main():
    print("--- Influenza Genome Motif Scanner ---")

    # 1. Build Matrix
    ll_matrix = build_log_likelihood_matrix(motifs)
    print("1. PWM constructed.")

    # 2. Find FASTA files in the 'genomes' subfolder
    search_path = os.path.join(GENOME_FOLDER, "*.fasta")
    fasta_files = glob.glob(search_path)

    if not fasta_files:
        print(f"No .fasta files found in the '{GENOME_FOLDER}' directory.")
        print(f"Current working directory: {os.getcwd()}")
        return

    print(f"2. Found {len(fasta_files)} genomes in '{GENOME_FOLDER}/'.")

    # 3. Process
    for file_path in fasta_files:
        print(f"\nProcessing: {os.path.basename(file_path)}...")

        genome_seq = read_fasta_file(file_path)
        if not genome_seq:
            continue

        if len(genome_seq) < MOTIF_LENGTH:
            print(f"   Skipping (Too short)")
            continue

        scores = scan_genome(genome_seq, ll_matrix)

        max_score = max(scores)
        max_index = scores.index(max_score)
        print(f"   Highest Signal: {max_score:.2f} at index {max_index}")

        plot_genome_signal(file_path, scores)

    print("\n--- Done. Charts have been saved in the main directory. ---")


if __name__ == "__main__":
    main()