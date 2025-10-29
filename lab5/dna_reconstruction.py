import random
import requests
from Bio import Entrez, SeqIO
from collections import defaultdict
import numpy as np

# Set your email for NCBI (required)
Entrez.email = "your.email@example.com"


def fetch_dna_sequence_from_ncbi(accession_id="NC_000913.3", start=1000, length=2000):
    """
    Fetch a DNA sequence from NCBI within the specified length range.
    Default uses E. coli genome as an example.
    """
    try:
        # Fetch sequence from NCBI
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta",
                               retmode="text", seq_start=start, seq_stop=start + length - 1)
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return str(record.seq)
    except:
        # Fallback: generate a random DNA sequence for demonstration
        print("NCBI fetch failed, generating random sequence for demonstration")
        bases = ['A', 'T', 'G', 'C']
        return ''.join(random.choices(bases, k=length))


def generate_random_samples(dna_sequence, num_samples=2000, sample_length_range=(100, 150)):
    """
    Generate random samples from the DNA sequence.
    """
    samples = []
    seq_length = len(dna_sequence)

    for _ in range(num_samples):
        # Random sample length between 100-150
        sample_length = random.randint(*sample_length_range)

        # Ensure we don't go beyond sequence boundaries
        max_start = max(0, seq_length - sample_length)
        start_pos = random.randint(0, max_start)

        sample = dna_sequence[start_pos:start_pos + sample_length]
        samples.append({
            'sequence': sample,
            'start_pos': start_pos,
            'end_pos': start_pos + sample_length
        })

    return samples


def reconstruct_sequence_overlap_method(samples, original_length):
    """
    Attempt to reconstruct the original sequence using overlap detection.
    This is a simplified approach to demonstrate the concept.
    """
    # Sort samples by their start position (cheating for demonstration)
    # In reality, we wouldn't know the positions
    samples_with_pos = [(s['start_pos'], s['sequence']) for s in samples]
    samples_with_pos.sort()

    reconstructed = ""
    last_end = 0

    for start_pos, sequence in samples_with_pos:
        if start_pos <= last_end:
            # There's an overlap, try to merge
            overlap_start = last_end - start_pos
            if overlap_start < len(sequence):
                reconstructed += sequence[overlap_start:]
                last_end = start_pos + len(sequence)
        else:
            # Gap in coverage
            reconstructed += 'N' * (start_pos - last_end)  # Fill gap with N's
            reconstructed += sequence
            last_end = start_pos + len(sequence)

    return reconstructed


def find_overlaps(samples, min_overlap=10):
    """
    Find overlaps between samples (without knowing positions).
    This is computationally expensive and demonstrates the main problem.
    """
    overlaps = []

    for i, sample1 in enumerate(samples):
        for j, sample2 in enumerate(samples[i + 1:], i + 1):
            seq1 = sample1['sequence']
            seq2 = sample2['sequence']

            # Check if seq1 overlaps with seq2
            for k in range(len(seq1) - min_overlap + 1):
                suffix = seq1[k:]
                if seq2.startswith(suffix) and len(suffix) >= min_overlap:
                    overlaps.append({
                        'sample1_idx': i,
                        'sample2_idx': j,
                        'overlap_length': len(suffix),
                        'type': 'seq1_to_seq2'
                    })

            # Check if seq2 overlaps with seq1
            for k in range(len(seq2) - min_overlap + 1):
                suffix = seq2[k:]
                if seq1.startswith(suffix) and len(suffix) >= min_overlap:
                    overlaps.append({
                        'sample1_idx': j,
                        'sample2_idx': i,
                        'overlap_length': len(suffix),
                        'type': 'seq2_to_seq1'
                    })

    return overlaps


def analyze_coverage(samples, original_length):
    """
    Analyze how well the samples cover the original sequence.
    """
    coverage = [0] * original_length

    for sample in samples:
        start = sample['start_pos']
        end = sample['end_pos']
        for i in range(start, min(end, original_length)):
            coverage[i] += 1

    coverage_stats = {
        'positions_covered': sum(1 for x in coverage if x > 0),
        'total_positions': original_length,
        'coverage_percentage': (sum(1 for x in coverage if x > 0) / original_length) * 100,
        'average_depth': np.mean(coverage),
        'positions_with_no_coverage': sum(1 for x in coverage if x == 0)
    }

    return coverage_stats, coverage


def main():
    print("DNA Sequence Reconstruction Demonstration")
    print("=" * 50)

    # Step 1: Get DNA sequence
    print("1. Fetching DNA sequence...")
    original_sequence = fetch_dna_sequence_from_ncbi()
    print(f"Original sequence length: {len(original_sequence)}")
    print(f"\nORIGINAL SEQUENCE:")
    print(f"First 200 bases: {original_sequence[:200]}")
    print(f"Last 200 bases:  {original_sequence[-200:]}")
    print(f"Full sequence: {original_sequence}")

    # Step 2: Generate random samples
    print("\n" + "=" * 50)
    print("2. Generating random samples...")
    samples = generate_random_samples(original_sequence)
    print(f"Generated {len(samples)} samples")
    print(
        f"Sample lengths range: {min(len(s['sequence']) for s in samples)} - {max(len(s['sequence']) for s in samples)}")

    # Show first 10 samples with their positions
    print(f"\nFIRST 10 SAMPLES:")
    for i, sample in enumerate(samples[:10]):
        print(
            f"Sample {i + 1:2d}: pos {sample['start_pos']:4d}-{sample['end_pos']:4d} (len={len(sample['sequence']):3d})")
        print(f"           {sample['sequence'][:50]}{'...' if len(sample['sequence']) > 50 else ''}")

    # Step 3: Analyze coverage
    print("\n" + "=" * 50)
    print("3. Analyzing coverage...")
    coverage_stats, coverage = analyze_coverage(samples, len(original_sequence))
    print(f"Coverage statistics:")
    for key, value in coverage_stats.items():
        print(f"  {key}: {value}")

    # Show coverage gaps
    uncovered_positions = [i for i, cov in enumerate(coverage) if cov == 0]
    if uncovered_positions:
        print(f"\nUNCOVERED POSITIONS: {uncovered_positions}")
        for pos in uncovered_positions:
            context_start = max(0, pos - 10)
            context_end = min(len(original_sequence), pos + 10)
            print(f"   Position {pos}: ...{original_sequence[context_start:context_end]}...")
    else:
        print("\nAll positions are covered!")

    # Show coverage depth distribution
    coverage_counts = {}
    for depth in coverage:
        coverage_counts[depth] = coverage_counts.get(depth, 0) + 1
    print(f"\nCOVERAGE DEPTH DISTRIBUTION:")
    for depth in sorted(coverage_counts.keys())[:10]:  # Show first 10 depths
        print(f"   Depth {depth:3d}: {coverage_counts[depth]:4d} positions")

    # Step 4: Attempt reconstruction (simplified)
    print("\n" + "=" * 50)
    print("4. RECONSTRUCTION RESULTS:")

    reconstructed = reconstruct_sequence_overlap_method(samples, len(original_sequence))

    print(f"Original length:      {len(original_sequence)}")
    print(f"Reconstructed length: {len(reconstructed)}")

    # Compare sequences
    if len(reconstructed) == len(original_sequence):
        matches = sum(1 for i, (a, b) in enumerate(zip(original_sequence, reconstructed)) if a == b)
        accuracy = (matches / len(original_sequence)) * 100
        print(f"Reconstruction accuracy: {accuracy:.2f}%")

        # Show first mismatches
        mismatches = [(i, a, b) for i, (a, b) in enumerate(zip(original_sequence, reconstructed)) if a != b]
        if mismatches:
            print(f"\nFIRST 5 MISMATCHES:")
            for i, (pos, orig, recon) in enumerate(mismatches[:5]):
                print(f"   Position {pos}: Original='{orig}' vs Reconstructed='{recon}'")
    else:
        print(f"Length mismatch - Original: {len(original_sequence)}, Reconstructed: {len(reconstructed)}")

    # Show comparison of sequences
    print(f"\nSEQUENCE COMPARISON (first 200 bases):")
    print(f"Original     : {original_sequence[:200]}")
    print(f"Reconstructed: {reconstructed[:200] if len(reconstructed) >= 200 else reconstructed}")

    if len(original_sequence) > 200:
        print(f"\nSEQUENCE COMPARISON (last 200 bases):")
        print(f"Original     : {original_sequence[-200:]}")
        print(f"Reconstructed: {reconstructed[-200:] if len(reconstructed) >= 200 else 'TOO SHORT'}")

    # Demonstrate the overlap finding problem
    print("\n" + "=" * 50)
    print("5. Demonstrating overlap detection complexity...")
    small_sample = samples[:50]  # Use only 50 samples for demo
    overlaps = find_overlaps(small_sample)
    print(f"Found {len(overlaps)} potential overlaps in {len(small_sample)} samples")
    print("This is O(n²) complexity and gets very expensive with 2000 samples!")

    # Show first few overlaps
    if overlaps:
        print(f"\nFIRST 5 OVERLAPS DETECTED:")
        for i, overlap in enumerate(overlaps[:5]):
            sample1_idx = overlap['sample1_idx']
            sample2_idx = overlap['sample2_idx']
            sample1 = small_sample[sample1_idx]
            sample2 = small_sample[sample2_idx]
            print(f"\nOverlap {i + 1}:")
            print(f"  Sample {sample1_idx} (pos {sample1['start_pos']}-{sample1['end_pos']})")
            print(f"  Sample {sample2_idx} (pos {sample2['start_pos']}-{sample2['end_pos']})")
            print(f"  Overlap length: {overlap['overlap_length']} bases")
            print(f"  Type: {overlap['type']}")

    # Calculate computational complexity
    full_comparisons = (len(samples) * (len(samples) - 1)) // 2
    subset_comparisons = (len(small_sample) * (len(small_sample) - 1)) // 2
    scaling_factor = full_comparisons / subset_comparisons

    print(f"\nCOMPUTATIONAL COMPLEXITY ANALYSIS:")
    print(f"  Subset ({len(small_sample)} samples): {subset_comparisons:,} comparisons")
    print(f"  Full set ({len(samples)} samples): {full_comparisons:,} comparisons")
    print(f"  Scaling factor: {scaling_factor:.0f}x more comparisons needed!")
    print(f"  Estimated overlaps for full set: {len(overlaps) * scaling_factor:.0f}")

    print(f"\nKEY PROBLEMS DEMONSTRATED:")
    print(f"  1. Length mismatch despite good coverage")
    print(f"  2. O(n²) computational complexity")
    print(f"  3. Ambiguous overlaps (many false positives)")
    print(f"  4. No position information in real scenarios")
    print(f"  5. Error propagation in reconstruction")


if __name__ == "__main__":
    main()