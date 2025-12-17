import math

motifs = [
    "GAGGTAAAC",
    "TCCGTAAGT",
    "CAGGTTGGA",
    "ACAGTCAGT",
    "TAGGTCATT",
    "TAGGTACTG",
    "ATGGTAACT",
    "CAGGTATAC",
    "TGTGTGAGT",
    "AAGGTAAGT"
]

sequence_s = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"

motif_length = 9
num_sequences = len(motifs)
bases = ['A', 'C', 'G', 'T']
background_freq = 0.25


count_matrix = {b: [0]*motif_length for b in bases}
prob_matrix = {b: [0.0]*motif_length for b in bases}
log_likelihood_matrix = {b: [0.0]*motif_length for b in bases}

for seq in motifs:
    for i, char in enumerate(seq):
        count_matrix[char][i] += 1

for b in bases:
    for i in range(motif_length):
        prob_matrix[b][i] = count_matrix[b][i] / num_sequences


for b in bases:
    for i in range(motif_length):
        p_obs = prob_matrix[b][i]
        if p_obs > 0:
            score = math.log(p_obs / background_freq)
        else:
            score = -99.0  

def print_matrix(title, matrix):
    print(f"\n--- {title} ---")
    print("   " + " ".join(f"{i+1:5}" for i in range(motif_length)))
    for b in bases:
        row_vals = " ".join(f"{val:5}" for val in matrix[b])
        print(f"{b} | {row_vals}")

print_matrix("Count Matrix", count_matrix)
print_matrix("Probability Matrix", prob_matrix)
print_matrix("Log-Likelihood Matrix", log_likelihood_matrix)


print(f"\n--- Analyzing Sequence S ---")
print(f"Sequence: {sequence_s}\n")

best_score = -999
best_window = ""
best_index = -1

print(f"{'Index':<6} {'Window':<12} {'Score':<6}")
print("-" * 30)

for i in range(len(sequence_s) - motif_length + 1):
    window = sequence_s[i : i + motif_length]
    current_score = 0
    
    for pos, char in enumerate(window):
        current_score += log_likelihood_matrix[char][pos]
        
    if current_score > -20:  
        print(f"{i:<6} {window:<12} {current_score:.2f}")
        
    if current_score > best_score:
        best_score = current_score
        best_window = window
        best_index = i

print("\n--- Conclusion ---")
print(f"Best match found at index {best_index}: {best_window}")
print(f"Score: {best_score:.2f}")

if best_score > 0:
    print("Result: SIGNAL DETECTED. This is likely an exon-intron border.")
else:
    print("Result: No significant signal found.")