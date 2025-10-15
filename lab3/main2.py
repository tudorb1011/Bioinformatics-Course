import os
import math
import matplotlib.pyplot as plt

def parse_fasta(file_path: str) -> str:

    if not os.path.exists(file_path):
        print(f"file not found at '{file_path}'")
        return ""

    sequence = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            sequence.append(line.strip())

    return "".join(sequence).upper()


def calculate_tm_values(dna_window: str, na_concentration: float):

    length = len(dna_window)
    g_count = dna_window.count('G')
    c_count = dna_window.count('C')
    a_count = dna_window.count('A')
    t_count = dna_window.count('T')
    gc_count = g_count + c_count

    # Formula 1: Basic Tm
    tm_basic = (4 * gc_count) + (2 * (a_count + t_count))

    # Formula 2: Salt-Adjusted Tm
    if length > 0:
        gc_percent = (gc_count / length) * 100
        tm_salt_adjusted = 81.5 + (16.6 * math.log10(na_concentration)) + (0.41 * gc_percent) - (600 / length)
    else:
        tm_salt_adjusted = 0

    return tm_basic, tm_salt_adjusted


def sliding_window_analysis(sequence: str, window_size: int, na_concentration: float):

    results = []
    if len(sequence) < window_size:
        print("seq is shorter than the window size")
        return results

    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        tm_basic, tm_salt = calculate_tm_values(window, na_concentration)
        results.append({
            "position": i + 1,
            "tm_basic": tm_basic,
            "tm_salt_adjusted": tm_salt
        })
    return results


def plot_tm_results(results, filename="tm_plot.png"):

    if not results:
        print("No data to plot.")
        return

    positions = [r['position'] for r in results]
    tm_basic_values = [r['tm_basic'] for r in results]
    tm_salt_values = [r['tm_salt_adjusted'] for r in results]

    plt.style.use('seaborn-v0_8-whitegrid')
    plt.figure(figsize=(16, 8))

    plt.plot(positions, tm_basic_values, label="Basic Formula: 4(G+C) + 2(A+T)", marker='o', linestyle='-',
             color='royalblue')
    plt.plot(positions, tm_salt_values, label="Salt-Adjusted Formula", marker='x', linestyle='--', color='crimson')

    plt.title("Sliding Window Melting Temperature (Tm) Analysis", fontsize=16)
    plt.xlabel("Window Start Position (bp)", fontsize=12)
    plt.ylabel("Melting Temperature (°C)", fontsize=12)
    plt.legend(fontsize=10)
    plt.grid(True)

    plt.xticks(range(min(positions), max(positions) + 1, max(1, len(positions) // 20)))
    plt.tight_layout()

    plt.savefig(filename)
    print(f"\nPlot successfully saved to '{filename}'")



if __name__ == "__main__":

    #define parameters
    fasta_file = "myfasta.fasta"
    NA_CONCENTRATION = 0.1  # Molar (e.g., 100 mM)
    WINDOW_SIZE = 8

    #read file
    dna_sequence = parse_fasta(fasta_file)

    #perform analysis if sequence is valid
    if dna_sequence:
        analysis_data = sliding_window_analysis(dna_sequence, WINDOW_SIZE, NA_CONCENTRATION)
        plot_tm_results(analysis_data)