import itertools
import tkinter as tk
from tkinter import ttk
from tkinter.scrolledtext import ScrolledText

DEFAULT_SEQUENCE = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"
ALPHABET = "ACGT"

def generate_kmers(k: int):
    return ["".join(p) for p in itertools.product(ALPHABET, repeat=k)]

def count_overlapping(seq: str, kmer: str) -> int:
    k = len(kmer)
    n = len(seq)
    if k > n:
        return 0
    return sum(1 for i in range(n - k + 1) if seq[i:i + k] == kmer)

def compute_percentages(seq: str, k: int):
    seq = seq.upper()
    total_windows = max(0, len(seq) - k + 1)
    kmers = generate_kmers(k)
    results = []
    for kmer in kmers:
        c = count_overlapping(seq, kmer)
        pct = (c / total_windows * 100.0) if total_windows else 0.0
        results.append((kmer, pct))
    return sorted(results, key=lambda x: x[0])

def format_results(title: str, data):
    lines = [title]
    for kmer, pct in data:
        lines.append(f"{kmer}: {pct:.2f}%")
    return "\n".join(lines)

def run_analysis(seq: str) -> str:
    di = compute_percentages(seq, 2)
    tri = compute_percentages(seq, 3)
    di_text = format_results("Dinucleotide percentages (k=2):", di)
    tri_text = format_results("Trinucleotide percentages (k=3):", tri)
    return f"{di_text}\n\n{tri_text}\n"

def build_gui():
    root = tk.Tk()
    root.title("K-mer Percentage Analyzer")

    main = ttk.Frame(root, padding=12)
    main.grid(sticky="nsew")
    root.rowconfigure(0, weight=1)
    root.columnconfigure(0, weight=1)
    main.rowconfigure(3, weight=1)
    main.columnconfigure(0, weight=1)


    seq_label = ttk.Label(main, text="DNA sequence:")
    seq_label.grid(row=0, column=0, sticky="w")
    seq_var = tk.StringVar(value=DEFAULT_SEQUENCE)
    seq_entry = ttk.Entry(main, textvariable=seq_var, width=80)
    seq_entry.grid(row=1, column=0, sticky="ew", pady=(0, 8))


    out = ScrolledText(main, width=90, height=24, wrap="none")
    out.configure(font=("Courier New", 10))
    out.grid(row=3, column=0, sticky="nsew", pady=(8, 0))


    def on_run(event=None):
        btn.config(state="disabled")
        out.delete("1.0", "end")
        seq = seq_var.get().strip()
        if not seq:
            out.insert("end", "Please enter a non-empty DNA sequence.\n")
            btn.config(state="normal")
            return
        result = run_analysis(seq)
        out.insert("end", result)
        btn.config(state="normal")

    btn = ttk.Button(main, text="Run analysis", command=on_run)
    btn.grid(row=2, column=0, sticky="w")

    root.bind("<Return>", on_run)
    seq_entry.focus_set()
    return root

if __name__ == "__main__":
    app = build_gui()
    app.mainloop()
