import random
import urllib.request
import time
import os

def get_reverse_complement(seq):
    """Returns the reverse complement of a DNA sequence."""
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(complement_map.get(base, 'N') for base in reversed(seq))

def generate_random_dna(length):
    """Generates a random DNA sequence of given length."""
    return "".join(random.choice("ATCG") for _ in range(length))

class TransposonSimulation:
    def __init__(self, target_len_min=200, target_len_max=400):
        self.seq_length = random.randint(target_len_min, target_len_max)
        self.sequence = generate_random_dna(self.seq_length)
        self.transposons = [] # Stores (start, end, ir_seq)

    def insert_transposons(self, count=3):

        print(f"\n--- Generating Artificial Sequence ({self.seq_length} bp) ---")
        
        for i in range(count):
            ir_len = random.randint(4, 6) 
            payload_len = random.randint(10, 30) 
            
            ir_seq = generate_random_dna(ir_len)
            ir_rc = get_reverse_complement(ir_seq)
            payload = generate_random_dna(payload_len)
            
            full_transposon = ir_seq + payload + ir_rc
            
            # 2. Insert into random position
            if len(self.sequence) == 0:
                insert_pos = 0
            else:
                insert_pos = random.randint(0, len(self.sequence))
            
            self.sequence = (
                self.sequence[:insert_pos] + 
                full_transposon + 
                self.sequence[insert_pos:]
            )
            
            print(f"Inserted TE #{i+1} at index {insert_pos}: IR='{ir_seq}' Len={len(full_transposon)}")
            
        self.seq_length = len(self.sequence)
        print(f"Final Sequence Length: {self.seq_length}")
        print(f"Sequence: {self.sequence}")


def find_transposons(dna_sequence, min_ir=4, max_ir=6, min_dist=10, max_dist=100):
    detected = []
    seq_len = len(dna_sequence)
    
    for i in range(seq_len):
        for k in range(min_ir, max_ir + 1):
            if i + k > seq_len: continue
            
            start_ir = dna_sequence[i : i+k]
            expected_end_ir = get_reverse_complement(start_ir)
            
            search_start = i + k + min_dist
            search_end = min(i + k + max_dist + k, seq_len)
            
            if search_start >= seq_len: continue
            
            window = dna_sequence[search_start : search_end]
            match_pos_in_window = window.find(expected_end_ir)
            
            if match_pos_in_window != -1:
                end_ir_start_index = search_start + match_pos_in_window
                end_ir_end_index = end_ir_start_index + k
                
                detected.append({
                    "start": i,
                    "end": end_ir_end_index,
                    "ir_seq": start_ir,
                    "ir_len": k,
                    "payload_len": end_ir_start_index - (i + k)
                })
                
    return detected


def get_genome_sequence(filename, accession):
    if os.path.exists(filename):
        print(f"Reading local file: {filename}")
        try:
            with open(filename, 'r') as f:
                data = f.read()
        except IOError as e:
            print(f"Error reading file {filename}: {e}")
            return ""
    else:
     
        try:
            with urllib.request.urlopen(url) as response:
                data = response.read().decode('utf-8')
                with open(filename, 'w') as f:
                    f.write(data)
                print(f"Saved data to {filename}")
        except Exception as e:
            print(f"Error downloading {accession}: {e}")
            return ""

    lines = data.split('\n')
    sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
    return sequence

def analyze_real_genomes():
    target_files = [
        ("C:\\Users\\RAUL\\OneDrive\\Desktop\\Politehnica\\Anul4\\Bioinformatics\\Project_L8\\L8\\Escherichia.fasta", "NC_000913"),     
        ("C:\\Users\\RAUL\\OneDrive\\Desktop\\Politehnica\\Anul4\\Bioinformatics\\Project_L8\\L8\\Bacillus.fasta", "NC_000964"),        
        ("C:\\Users\\RAUL\\OneDrive\\Desktop\\Politehnica\\Anul4\\Bioinformatics\\Project_L8\\L8\\Staphylococcus.fasta", "NC_007795")  
    ]
    
    print("\n" + "="*40)
    print("PART 3: REAL GENOME ANALYSIS (FROM FILES)")
    print("="*40)
    
    for filename, accession in target_files:
        print(f"\nProcessing: {filename}")
        
        seq = get_genome_sequence(filename, accession)
        if not seq: 
            print("Skipping (No sequence found).")
            continue
            
        print(f"Total Sequence Length: {len(seq)} bp")
        
        scan_limit = 5000 
        sub_seq = seq[:scan_limit]
        
        print(f"Scanning first {scan_limit} bp for possible transposons...")
        print("(Criteria: IR len 4-6, Payload 50-1000 bp)")
        
        hits = find_transposons(sub_seq, min_ir=4, max_ir=6, min_dist=50, max_dist=1000)
        
        print(f"Found {len(hits)} potential inverted repeat pairs.")
        if hits:
            print("Top 5 examples:")
            for i, hit in enumerate(hits[:5]):
                print(f"  {i+1}. Loc: {hit['start']}-{hit['end']} | IR: {hit['ir_seq']} | Payload: {hit['payload_len']}bp")
        else:
            print("No hits found in the scanned region.")
            
        time.sleep(1)


if __name__ == "__main__":
    print("="*40)
    print("PART 1 & 2: ARTIFICIAL SIMULATION")
    print("="*40)
    
    sim = TransposonSimulation()
    sim.insert_transposons(count=3)
    
    print("\n--- Running Detection Software ---")
    results = find_transposons(sim.sequence, min_ir=4, max_ir=6, min_dist=10, max_dist=100)
    
    print(f"Detected {len(results)} potential elements:")
    for res in results:
        print(f"  Found at {res['start']}-{res['end']} (IR: {res['ir_seq']}, Payload: {res['payload_len']}bp)")

    analyze_real_genomes()