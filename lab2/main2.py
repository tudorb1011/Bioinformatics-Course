from typing import List, Tuple, Dict

def scan_kmers(seq: str, k: int) -> List[Tuple[str, int, float]]:

    n = len(seq)
    if k <= 0:
        raise ValueError("k must be positive")
    total_windows = max(0, n - k + 1)
    if total_windows == 0:
        return []

    counts: Dict[str, int] = {}
    first_idx: Dict[str, int] = {}

    for i in range(total_windows):
        kmer = seq[i:i + k]
        if kmer not in counts:
            counts[kmer] = 0
            first_idx[kmer] = i
        counts[kmer] += 1

    ordered = sorted(counts.keys(), key=lambda t: first_idx[t])
    return [(kmer, counts[kmer], counts[kmer] / total_windows * 100.0) for kmer in ordered]

def print_results(seq: str) -> None:
    print(f"Sequence: {seq}")
    di = scan_kmers(seq, 2)
    tri = scan_kmers(seq, 3)

    print("\nDinucleotides present (in order of first appearance):")
    for kmer, cnt, pct in di:
        print(f"{kmer}: count={cnt}, pct={pct:.2f}%")

    print("\nTrinucleotides present (in order of first appearance):")
    for kmer, cnt, pct in tri:
        print(f"{kmer}: count={cnt}, pct={pct:.2f}%")

if __name__ == "__main__":
    # Example: s = 'abaa'
    print_results("abaa")
    # You can also test with a DNA sequence:
    # print_results("TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA")
