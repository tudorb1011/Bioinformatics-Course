import math

def calculate_tm(dna_sequence: str, na_concentration: float = 0.1):

    sequence = dna_sequence.upper()
    valid_bases = 'ATGC'

    if not all(base in valid_bases for base in sequence):
        return {"error": "Please use only A, T, G, C"}

    length = len(sequence)
    if length == 0:
        return {"error": "DNA sequence cant be empty"}


    a_count = sequence.count('A')
    t_count = sequence.count('T')
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    gc_count = g_count + c_count



    # Formula 1 = Basic Calculation
    # Tm = 4(G + C) + 2(A + T) °C

    tm_basic = (4 * gc_count) + (2 * (a_count + t_count))

    # Formula 2 = Salt-Adjusted Calculation
    # Tm = 81.5 + 16.6(log10([Na+])) + 0.41*(%GC) – 600/length

    gc_percent = (gc_count / length) * 100
    tm_salt_adjusted = 81.5 + (16.6 * math.log10(na_concentration)) + (0.41 * gc_percent) - (600 / length)


    return {
        "sequence": sequence,
        "length": length,
        "gc_content_%": f"{gc_percent:.2f}",
        "tm_basic_celsius": f"{tm_basic:.2f}",
        "tm_salt_adjusted_celsius": f"{tm_salt_adjusted:.2f}"
    }


my_dna_primer = "AGCTTAGGCTAC"
primer_na_concentration = 0.15  # Molar (M)

results = calculate_tm(my_dna_primer, primer_na_concentration)


if "error" in results:
    print(f"Error: {results['error']}")
else:
    print(f"--- Tm result ---")
    print(f"Sequence:           {results['sequence']}")
    print(f"Length:             {results['length']}")
    print(f"GC Content:         {results['gc_content_%']}%")
    print("-" * 30)
    print(f"Basic Tm:           {results['tm_basic_celsius']} °C")
    print(f"Salt-Adjusted Tm:   {results['tm_salt_adjusted_celsius']} °C (at {primer_na_concentration}M Na+)")
    print("-" * 30)