import sys
import os
import glob
import matplotlib.pyplot as plt

def get_filename():
    if len(sys.argv) > 1:
        return sys.argv[1]
    # Auto-pick the first FASTA from the current folder
    for pattern in ("*.fasta", "*.fa"):
        matches = sorted(glob.glob(pattern))
        if matches:
            return matches[0]
    print("No FASTA file found in the current folder and none provided as an argument.")
    sys.exit(1)

filename = get_filename()
window_size = 30
smooth_window = 10

all_percentages = []

with open(filename, 'r') as f:
    for line in f:
        trimmed = line.strip()
        if trimmed.startswith('>'):
            continue

        for i in range(len(trimmed) - window_size + 1):
            chunk = trimmed[i:i + window_size]
            occurrences = {}
            for letter in chunk:
                occurrences[letter] = occurrences.get(letter, 0) + 1

            percentages = {k: v / len(chunk) for k, v in occurrences.items()}
            all_percentages.append(percentages)

bases = ['A', 'C', 'G', 'T']
window_indices = list(range(1, len(all_percentages) + 1))

freq_data = {base: [] for base in bases}
for pct in all_percentages:
    for base in bases:
        freq_data[base].append(pct.get(base, 0))

def smooth(data, size):
    smoothed = []
    for i in range(len(data)):
        start = max(0, i - size + 1)
        smoothed.append(sum(data[start:i+1]) / (i - start + 1))
    return smoothed

plt.figure(figsize=(12, 6))
for base in bases:
    smoothed_data = smooth(freq_data[base], smooth_window)
    plt.plot(window_indices, smoothed_data, label=base)

plt.xlabel("Window number")
plt.ylabel("Frequency")
plt.title(f"Smoothed Base Frequencies per {window_size}-base Sliding Window")
plt.legend()
plt.grid(True)
plt.show()