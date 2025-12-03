import re
import sys
import os

def read_fasta(filename):

    if not os.path.exists(filename):
        print(f"Error: The file '{filename}' was not found.")
        print("Please ensure 'sequence.fasta' is in the same directory as this script.")
        sys.exit(1)

    with open(filename, 'r') as f:
        lines = f.readlines()

    sequence = ""
    for line in lines:
        line = line.strip()
        if not line.startswith(">"):
            sequence += line
            
    return sequence.upper() 

enzymes = {
    "EcoRI":   {"seq": "GAATTC", "offset": 1},
    "BamHI":   {"seq": "GGATCC", "offset": 1},
    "HindIII": {"seq": "AAGCTT", "offset": 1},
    "TaqI":    {"seq": "TCGA",   "offset": 1},
    "HaeIII":  {"seq": "GGCC",   "offset": 2}
}

def digest_dna(dna, enzyme_info, is_circular=True):
    seq_pattern = enzyme_info["seq"]
    cut_offset = enzyme_info["offset"]
    
    matches = [m.start() for m in re.finditer(f'(?={seq_pattern})', dna)]
    
    cut_sites = [m + cut_offset for m in matches]
    cut_sites.sort()
    
    num_cleavages = len(cut_sites)
    fragments = []
    
    if num_cleavages == 0:
        fragments = [len(dna)] 
    else:
        if is_circular:
            for i in range(len(cut_sites) - 1):
                fragments.append(cut_sites[i+1] - cut_sites[i])
            wrap_fragment = (len(dna) - cut_sites[-1]) + cut_sites[0]
            fragments.append(wrap_fragment)
        else:
            current_pos = 0
            for cut in cut_sites:
                fragments.append(cut - current_pos)
                current_pos = cut
            fragments.append(len(dna) - current_pos)
            
    fragments.sort(reverse=True) 
    return num_cleavages, cut_sites, fragments

def simulate_gel(results):
    print("\n" + "="*65)
    print("      SIMULATED ELECTROPHORESIS GEL (1.5% Agarose)")
    print("="*65)
    
    resolution_bins = [
        3000, 2500, 2000, 1500, 1200, 1000, 
        900, 800, 700, 600, 500, 400, 300, 200, 100, 50
    ]
    
    header = "Size(bp)|"
    for name in results.keys():
        header += f" {name.center(7)} |"
    print(header)
    print("-" * len(header))
    
    for limit in resolution_bins:
        row_str = f"{str(limit).rjust(7)} |"
        
        for name, data in results.items():
            fragments = data['fragments']
            
            in_bin = []
            for f in fragments:
                if limit >= f > (limit * 0.75): 
                     in_bin.append(f)
            
            if in_bin:
                if len(in_bin) == 1:
                    row_str += f"  [{in_bin[0]}]  |" 
                else:
                    row_str += "  [====] |" 
            else:
                row_str += "         |"
                
        print(row_str)
    print("-" * len(header))
    print("* Fragments < 50bp often run off the gel or are invisible.")


def main():
    filename = ("sequence.fasta")
    
    print(f"Reading '{filename}'...")
    dna_sequence = read_fasta(filename)
    
    print(f"Read {len(dna_sequence)} base pairs.")
    print(f"Topology Assumption: Circular (Plasmid)")
    print("-" * 40)

    all_results = {}

    for name, info in enzymes.items():
        cuts, sites, frags = digest_dna(dna_sequence, info, is_circular=True)
        all_results[name] = {'fragments': frags}
        
        print(f"\nEnzyme: {name}")
        print(f"  Sequence:  5'-{info['seq']}-3'")
        print(f"  Cleavages: {cuts}")
        if cuts > 0:
            print(f"  Positions: {sites}")
            print(f"  Fragments: {frags}")
        else:
            print("  Result:    No cuts (Supercoiled/Nicked DNA)")

    simulate_gel(all_results)

if __name__ == "__main__":
    main()