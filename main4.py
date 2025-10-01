import tkinter as tk
from tkinter import filedialog, messagebox
import os

def read_fasta(file_path: str) -> str:
    """Read a FASTA file and return the sequence content as a single string."""
    sequence = []
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            sequence.append(line.upper())
    return "".join(sequence)


def calculate_percentage(sequence: str) -> dict:
    """Calculate frequency and percentage of each symbol in the sequence."""
    alphabet = sorted(set(sequence))
    total_length = len(sequence)
    percentage_dict = {}

    for letter in alphabet:
        count = sequence.count(letter)
        percentage = (count / total_length) * 100
        percentage_dict[letter] = {
            'count': count,
            'percentage': round(percentage, 2)
        }

    return percentage_dict


# -------------------- GUI --------------------

def open_fasta_file():
    """Open a file dialog to select a FASTA file."""
    global fasta_file_path
    fasta_file_path = filedialog.askopenfilename(
        title="Select FASTA File",
        filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")]
    )
    if fasta_file_path:
        label_file.config(text=f"Selected file: {os.path.basename(fasta_file_path)}")
    else:
        label_file.config(text="No file selected.")


def analyze_sequence():
    """Read FASTA, calculate percentages, and display results."""
    if not fasta_file_path:
        messagebox.showerror("Error", "Please select a FASTA file first!")
        return

    try:
        seq = read_fasta(fasta_file_path)
        if not seq:
            messagebox.showerror("Error", "FASTA file is empty or invalid!")
            return

        result = calculate_percentage(seq)

        # Display results
        output_text.delete("1.0", tk.END)
        output_text.insert(tk.END, f"Alphabet: {', '.join(sorted(set(seq)))}\n")
        output_text.insert(tk.END, f"Total length: {len(seq)}\n\nComposition:\n")

        for letter, stats in result.items():
            output_text.insert(
                tk.END, f"Letter: {letter} | Count: {stats['count']} | Percentage: {stats['percentage']}%\n"
            )

    except Exception as e:
        messagebox.showerror("Error", f"Something went wrong:\n{e}")


# ---------- Build Window ----------
root = tk.Tk()
root.title("FASTA Sequence Analyzer")
root.geometry("600x500")
root.resizable(False, False)

fasta_file_path = ""

# Title label
title_label = tk.Label(root, text="FASTA Alphabet & Percentage Analyzer", font=("Arial", 16, "bold"))
title_label.pack(pady=10)

# Buttons
btn_open = tk.Button(root, text="📁 Open FASTA File", command=open_fasta_file, font=("Arial", 12))
btn_open.pack(pady=5)

label_file = tk.Label(root, text="No file selected.", font=("Arial", 11), fg="gray")
label_file.pack(pady=5)

btn_analyze = tk.Button(root, text="🔍 Analyze Sequence", command=analyze_sequence, font=("Arial", 12), bg="#4CAF50", fg="white")
btn_analyze.pack(pady=10)

# Output text area
output_text = tk.Text(root, height=20, width=70, font=("Courier", 11))
output_text.pack(pady=10)

root.mainloop()
