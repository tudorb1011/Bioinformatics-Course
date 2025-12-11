"""
Large Genome Local Alignment using Smith-Waterman Algorithm
with Sliding Window Approach for Influenza and COVID-19 Genomes

This implementation divides large genomes into manageable regions,
performs local alignment on each region, and visualizes the similarities.

Requirements: Python 3, tkinter (built-in), numpy
"""

import tkinter as tk
from tkinter import ttk, scrolledtext, filedialog, messagebox
import numpy as np
import os
from datetime import datetime

class GenomeLocalAlignment:
    def __init__(self, root):
        self.root = root
        self.root.title("Large Genome Local Alignment - Influenza vs COVID-19")
        self.root.geometry("1400x900")
        
        # Genome data
        self.genome1 = ""
        self.genome2 = ""
        self.genome1_name = ""
        self.genome2_name = ""
        
        # Alignment parameters
        self.match_score = 2
        self.mismatch_score = -1
        self.gap_penalty = -2
        self.window_size = 1000  # Size of each alignment region
        self.step_size = 500     # Overlap between windows
        
        # Results storage
        self.alignment_results = []
        
        self.create_ui()
        
    def create_ui(self):
        # Main container
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(4, weight=1)
        
        # Title
        title_label = tk.Label(main_frame, text="Large Genome Local Alignment System", 
                              font=("Arial", 16, "bold"))
        title_label.grid(row=0, column=0, columnspan=3, pady=10)
        
        # File loading section
        file_frame = ttk.LabelFrame(main_frame, text="Genome Files", padding="10")
        file_frame.grid(row=1, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)
        
        ttk.Button(file_frame, text="Load COVID-19 Genome", 
                  command=self.load_covid_genome).grid(row=0, column=0, padx=5)
        self.covid_label = ttk.Label(file_frame, text="No file loaded")
        self.covid_label.grid(row=0, column=1, padx=5, sticky=tk.W)
        
        ttk.Button(file_frame, text="Load Influenza Genome", 
                  command=self.load_influenza_genome).grid(row=1, column=0, padx=5, pady=5)
        self.influenza_label = ttk.Label(file_frame, text="No file loaded")
        self.influenza_label.grid(row=1, column=1, padx=5, sticky=tk.W)
        
        # Parameters section
        param_frame = ttk.LabelFrame(main_frame, text="Alignment Parameters", padding="10")
        param_frame.grid(row=2, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)
        
        ttk.Label(param_frame, text="Match Score:").grid(row=0, column=0, sticky=tk.W, padx=5)
        self.match_entry = ttk.Entry(param_frame, width=10)
        self.match_entry.insert(0, "2")
        self.match_entry.grid(row=0, column=1, padx=5)
        
        ttk.Label(param_frame, text="Mismatch Score:").grid(row=0, column=2, sticky=tk.W, padx=5)
        self.mismatch_entry = ttk.Entry(param_frame, width=10)
        self.mismatch_entry.insert(0, "-1")
        self.mismatch_entry.grid(row=0, column=3, padx=5)
        
        ttk.Label(param_frame, text="Gap Penalty:").grid(row=0, column=4, sticky=tk.W, padx=5)
        self.gap_entry = ttk.Entry(param_frame, width=10)
        self.gap_entry.insert(0, "-2")
        self.gap_entry.grid(row=0, column=5, padx=5)
        
        ttk.Label(param_frame, text="Window Size:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
        self.window_entry = ttk.Entry(param_frame, width=10)
        self.window_entry.insert(0, "1000")
        self.window_entry.grid(row=1, column=1, padx=5, pady=5)
        
        ttk.Label(param_frame, text="Step Size:").grid(row=1, column=2, sticky=tk.W, padx=5, pady=5)
        self.step_entry = ttk.Entry(param_frame, width=10)
        self.step_entry.insert(0, "500")
        self.step_entry.grid(row=1, column=3, padx=5, pady=5)
        
        # Control buttons
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=3, column=0, columnspan=3, pady=10)
        
        ttk.Button(button_frame, text="Start Alignment", 
                  command=self.start_alignment, 
                  style="Accent.TButton").pack(side=tk.LEFT, padx=5)
        
        ttk.Button(button_frame, text="Export Results", 
                  command=self.export_results).pack(side=tk.LEFT, padx=5)
        
        ttk.Button(button_frame, text="Clear Results", 
                  command=self.clear_results).pack(side=tk.LEFT, padx=5)
        
        # Results section with tabs
        notebook = ttk.Notebook(main_frame)
        notebook.grid(row=4, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)
        
        # Tab 1: Alignment results
        results_frame = ttk.Frame(notebook)
        notebook.add(results_frame, text="Alignment Results")
        
        self.results_text = scrolledtext.ScrolledText(results_frame, width=100, height=30, 
                                                      font=("Courier", 9))
        self.results_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Tab 2: Similarity visualization
        viz_frame = ttk.Frame(notebook)
        notebook.add(viz_frame, text="Similarity Visualization")
        
        self.viz_text = scrolledtext.ScrolledText(viz_frame, width=100, height=30, 
                                                  font=("Courier", 8))
        self.viz_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Tab 3: Statistics
        stats_frame = ttk.Frame(notebook)
        notebook.add(stats_frame, text="Statistics")
        
        self.stats_text = scrolledtext.ScrolledText(stats_frame, width=100, height=30, 
                                                    font=("Courier", 10))
        self.stats_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Status bar
        self.status_var = tk.StringVar()
        self.status_var.set("Ready")
        status_bar = ttk.Label(main_frame, textvariable=self.status_var, 
                              relief=tk.SUNKEN, anchor=tk.W)
        status_bar.grid(row=5, column=0, columnspan=3, sticky=(tk.W, tk.E))
        
    def load_fasta(self, filename):
        """Load a FASTA file and return the sequence"""
        try:
            with open(filename, 'r') as f:
                lines = f.readlines()
            
            sequence = ""
            header = ""
            
            for line in lines:
                line = line.strip()
                if line.startswith('>'):
                    header = line[1:]
                else:
                    sequence += line.upper()
            
            return sequence, header
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load file: {str(e)}")
            return "", ""
    
    def load_covid_genome(self):
        """Load COVID-19 genome file"""
        # First try to load covid19.fasta from current directory
        if os.path.exists("Covid19.fasta"):
            filename = "Covid19.fasta"
        else:
            filename = filedialog.askopenfilename(
                title="Select COVID-19 Genome File",
                filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")]
            )
        
        if filename:
            self.genome1, self.genome1_name = self.load_fasta(filename)
            if self.genome1:
                self.covid_label.config(text=f"Loaded: {os.path.basename(filename)} ({len(self.genome1)} bp)")
                self.status_var.set(f"COVID-19 genome loaded: {len(self.genome1)} base pairs")
    
    def load_influenza_genome(self):
        """Load Influenza genome file"""
        filename = filedialog.askopenfilename(
            title="Select Influenza Genome File",
            filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")]
        )
        
        if filename:
            self.genome2, self.genome2_name = self.load_fasta(filename)
            if self.genome2:
                self.influenza_label.config(text=f"Loaded: {os.path.basename(filename)} ({len(self.genome2)} bp)")
                self.status_var.set(f"Influenza genome loaded: {len(self.genome2)} base pairs")
    
    def smith_waterman_local(self, seq1, seq2, match, mismatch, gap):
        """
        Smith-Waterman algorithm for local sequence alignment
        Returns the scoring matrix and the maximum score with its position
        """
        n, m = len(seq1), len(seq2)
        
        # Initialize scoring matrix with zeros (local alignment)
        score_matrix = np.zeros((n + 1, m + 1), dtype=int)
        
        max_score = 0
        max_pos = (0, 0)
        
        # Fill the scoring matrix
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                # Calculate scores for three possible moves
                match_score = score_matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
                delete_score = score_matrix[i-1][j] + gap
                insert_score = score_matrix[i][j-1] + gap
                
                # Take maximum, but not less than 0 (local alignment)
                score_matrix[i][j] = max(0, match_score, delete_score, insert_score)
                
                # Track maximum score and position
                if score_matrix[i][j] > max_score:
                    max_score = score_matrix[i][j]
                    max_pos = (i, j)
        
        return score_matrix, max_score, max_pos
    
    def traceback_local(self, score_matrix, seq1, seq2, max_pos, match, mismatch, gap):
        """
        Traceback for Smith-Waterman local alignment
        Returns the aligned sequences
        """
        align1, align2 = "", ""
        i, j = max_pos
        
        # Traceback until we hit a zero
        while i > 0 and j > 0 and score_matrix[i][j] > 0:
            current_score = score_matrix[i][j]
            diagonal_score = score_matrix[i-1][j-1]
            up_score = score_matrix[i-1][j]
            left_score = score_matrix[i][j-1]
            
            # Check which direction we came from
            if current_score == diagonal_score + (match if seq1[i-1] == seq2[j-1] else mismatch):
                align1 = seq1[i-1] + align1
                align2 = seq2[j-1] + align2
                i -= 1
                j -= 1
            elif current_score == up_score + gap:
                align1 = seq1[i-1] + align1
                align2 = "-" + align2
                i -= 1
            else:
                align1 = "-" + align1
                align2 = seq2[j-1] + align2
                j -= 1
        
        return align1, align2, (i, j)
    
    def calculate_similarity(self, align1, align2):
        """Calculate similarity percentage between aligned sequences"""
        matches = sum(1 for a, b in zip(align1, align2) if a == b and a != '-')
        total = len(align1)
        return (matches / total * 100) if total > 0 else 0
    
    def sliding_window_alignment(self, genome1, genome2, window_size, step_size, 
                                 match, mismatch, gap):
        """
        Perform local alignment using sliding windows
        This allows us to handle large genomes by breaking them into regions
        """
        results = []
        
        # Calculate number of windows needed
        num_windows1 = max(1, (len(genome1) - window_size) // step_size + 1)
        num_windows2 = max(1, (len(genome2) - window_size) // step_size + 1)
        
        total_comparisons = num_windows1 * num_windows2
        current_comparison = 0
        
        self.status_var.set(f"Starting alignment: {total_comparisons} comparisons to perform...")
        self.root.update()
        
        # Iterate through windows of genome1
        for i in range(num_windows1):
            start1 = i * step_size
            end1 = min(start1 + window_size, len(genome1))
            window1 = genome1[start1:end1]
            
            # Iterate through windows of genome2
            for j in range(num_windows2):
                start2 = j * step_size
                end2 = min(start2 + window_size, len(genome2))
                window2 = genome2[start2:end2]
                
                current_comparison += 1
                
                # Update status every 10 comparisons
                if current_comparison % 10 == 0:
                    progress = (current_comparison / total_comparisons) * 100
                    self.status_var.set(f"Progress: {current_comparison}/{total_comparisons} ({progress:.1f}%)")
                    self.root.update()
                
                # Perform local alignment on these windows
                score_matrix, max_score, max_pos = self.smith_waterman_local(
                    window1, window2, match, mismatch, gap
                )
                
                # Only process if there's a significant alignment
                if max_score > 50:  # Threshold for significant alignment
                    align1, align2, start_pos = self.traceback_local(
                        score_matrix, window1, window2, max_pos, match, mismatch, gap
                    )
                    
                    similarity = self.calculate_similarity(align1, align2)
                    
                    result = {
                        'window1_start': start1,
                        'window1_end': end1,
                        'window2_start': start2,
                        'window2_end': end2,
                        'score': max_score,
                        'similarity': similarity,
                        'align1': align1,
                        'align2': align2,
                        'local_start1': start_pos[0],
                        'local_start2': start_pos[1]
                    }
                    
                    results.append(result)
        
        # Sort results by score (highest first)
        results.sort(key=lambda x: x['score'], reverse=True)
        
        return results
    
    def start_alignment(self):
        """Start the alignment process"""
        # Validate inputs
        if not self.genome1 or not self.genome2:
            messagebox.showerror("Error", "Please load both genome files first!")
            return
        
        try:
            # Get parameters
            self.match_score = int(self.match_entry.get())
            self.mismatch_score = int(self.mismatch_entry.get())
            self.gap_penalty = int(self.gap_entry.get())
            self.window_size = int(self.window_entry.get())
            self.step_size = int(self.step_entry.get())
        except ValueError:
            messagebox.showerror("Error", "Please enter valid numeric parameters!")
            return
        
        # Clear previous results
        self.clear_results()
        
        # Perform alignment
        self.status_var.set("Performing sliding window alignment...")
        self.root.update()
        
        self.alignment_results = self.sliding_window_alignment(
            self.genome1, self.genome2, 
            self.window_size, self.step_size,
            self.match_score, self.mismatch_score, self.gap_penalty
        )
        
        # Display results
        self.display_results()
        self.display_visualization()
        self.display_statistics()
        
        self.status_var.set(f"Alignment complete! Found {len(self.alignment_results)} significant regions")
    
    def display_results(self):
        """Display alignment results in the results tab"""
        self.results_text.delete(1.0, tk.END)
        
        self.results_text.insert(tk.END, "="*100 + "\n")
        self.results_text.insert(tk.END, "GENOME LOCAL ALIGNMENT RESULTS\n")
        self.results_text.insert(tk.END, "="*100 + "\n\n")
        
        self.results_text.insert(tk.END, f"Genome 1: {self.genome1_name[:80]}\n")
        self.results_text.insert(tk.END, f"Length: {len(self.genome1)} bp\n\n")
        
        self.results_text.insert(tk.END, f"Genome 2: {self.genome2_name[:80]}\n")
        self.results_text.insert(tk.END, f"Length: {len(self.genome2)} bp\n\n")
        
        self.results_text.insert(tk.END, f"Parameters:\n")
        self.results_text.insert(tk.END, f"  Match: {self.match_score}, Mismatch: {self.mismatch_score}, Gap: {self.gap_penalty}\n")
        self.results_text.insert(tk.END, f"  Window Size: {self.window_size}, Step Size: {self.step_size}\n\n")
        
        self.results_text.insert(tk.END, "="*100 + "\n\n")
        
        # Display top alignments
        num_display = min(20, len(self.alignment_results))
        
        for idx, result in enumerate(self.alignment_results[:num_display], 1):
            self.results_text.insert(tk.END, f"Alignment #{idx}\n")
            self.results_text.insert(tk.END, f"{'-'*100}\n")
            self.results_text.insert(tk.END, f"Region 1: {result['window1_start']}-{result['window1_end']} bp\n")
            self.results_text.insert(tk.END, f"Region 2: {result['window2_start']}-{result['window2_end']} bp\n")
            self.results_text.insert(tk.END, f"Score: {result['score']}\n")
            self.results_text.insert(tk.END, f"Similarity: {result['similarity']:.2f}%\n\n")
            
            # Display alignment in chunks
            align1 = result['align1']
            align2 = result['align2']
            
            # Show first 200 characters of alignment
            display_length = min(200, len(align1))
            
            self.results_text.insert(tk.END, f"Alignment (showing first {display_length} positions):\n")
            
            chunk_size = 60
            for i in range(0, display_length, chunk_size):
                chunk1 = align1[i:i+chunk_size]
                chunk2 = align2[i:i+chunk_size]
                
                # Create match line
                match_line = ""
                for a, b in zip(chunk1, chunk2):
                    if a == b and a != '-':
                        match_line += "|"
                    elif a != '-' and b != '-':
                        match_line += "."
                    else:
                        match_line += " "
                
                self.results_text.insert(tk.END, f"  {chunk1}\n")
                self.results_text.insert(tk.END, f"  {match_line}\n")
                self.results_text.insert(tk.END, f"  {chunk2}\n\n")
            
            if len(align1) > display_length:
                self.results_text.insert(tk.END, f"  ... (alignment continues for {len(align1)} total positions)\n\n")
            
            self.results_text.insert(tk.END, "\n")
    
    def display_visualization(self):
        """Display a visual representation of similarity across genomes"""
        self.viz_text.delete(1.0, tk.END)
        
        self.viz_text.insert(tk.END, "="*100 + "\n")
        self.viz_text.insert(tk.END, "SIMILARITY MAP\n")
        self.viz_text.insert(tk.END, "="*100 + "\n\n")
        
        # Create a grid representation
        grid_size = 50
        len1 = len(self.genome1)
        len2 = len(self.genome2)
        
        bin_size1 = max(1, len1 // grid_size)
        bin_size2 = max(1, len2 // grid_size)
        
        # Initialize grid with zeros
        grid = [[0 for _ in range(grid_size)] for _ in range(grid_size)]
        
        # Fill grid based on alignment results
        for result in self.alignment_results:
            bin1 = min(grid_size - 1, result['window1_start'] // bin_size1)
            bin2 = min(grid_size - 1, result['window2_start'] // bin_size2)
            
            # Weight by similarity score
            grid[bin1][bin2] = max(grid[bin1][bin2], result['similarity'])
        
        # Display grid
        self.viz_text.insert(tk.END, "Genome 2 →\n")
        self.viz_text.insert(tk.END, "↓ Genome 1\n\n")
        
        # Create legend
        symbols = {
            'high': '█',      # 70-100% similarity
            'medium': '▓',    # 40-70% similarity
            'low': '░',       # 10-40% similarity
            'none': '·'       # 0-10% similarity
        }
        
        self.viz_text.insert(tk.END, f"Legend: {symbols['high']} = High (70-100%)  ")
        self.viz_text.insert(tk.END, f"{symbols['medium']} = Medium (40-70%)  ")
        self.viz_text.insert(tk.END, f"{symbols['low']} = Low (10-40%)  ")
        self.viz_text.insert(tk.END, f"{symbols['none']} = None (0-10%)\n\n")
        
        # Draw grid
        for i in range(grid_size):
            row = ""
            for j in range(grid_size):
                similarity = grid[i][j]
                if similarity >= 70:
                    row += symbols['high']
                elif similarity >= 40:
                    row += symbols['medium']
                elif similarity >= 10:
                    row += symbols['low']
                else:
                    row += symbols['none']
            self.viz_text.insert(tk.END, row + "\n")
        
        self.viz_text.insert(tk.END, f"\nEach cell represents approximately {bin_size1} x {bin_size2} base pairs\n")
    
    def display_statistics(self):
        """Display statistical analysis of the alignment"""
        self.stats_text.delete(1.0, tk.END)
        
        self.stats_text.insert(tk.END, "="*80 + "\n")
        self.stats_text.insert(tk.END, "ALIGNMENT STATISTICS\n")
        self.stats_text.insert(tk.END, "="*80 + "\n\n")
        
        if not self.alignment_results:
            self.stats_text.insert(tk.END, "No significant alignments found.\n")
            return
        
        # Calculate statistics
        num_alignments = len(self.alignment_results)
        scores = [r['score'] for r in self.alignment_results]
        similarities = [r['similarity'] for r in self.alignment_results]
        
        avg_score = np.mean(scores)
        max_score = np.max(scores)
        min_score = np.min(scores)
        
        avg_similarity = np.mean(similarities)
        max_similarity = np.max(similarities)
        min_similarity = np.min(similarities)
        
        self.stats_text.insert(tk.END, f"Total Significant Alignments Found: {num_alignments}\n\n")
        
        self.stats_text.insert(tk.END, "Alignment Scores:\n")
        self.stats_text.insert(tk.END, f"  Maximum:  {max_score}\n")
        self.stats_text.insert(tk.END, f"  Minimum:  {min_score}\n")
        self.stats_text.insert(tk.END, f"  Average:  {avg_score:.2f}\n")
        self.stats_text.insert(tk.END, f"  Std Dev:  {np.std(scores):.2f}\n\n")
        
        self.stats_text.insert(tk.END, "Similarity Percentages:\n")
        self.stats_text.insert(tk.END, f"  Maximum:  {max_similarity:.2f}%\n")
        self.stats_text.insert(tk.END, f"  Minimum:  {min_similarity:.2f}%\n")
        self.stats_text.insert(tk.END, f"  Average:  {avg_similarity:.2f}%\n")
        self.stats_text.insert(tk.END, f"  Std Dev:  {np.std(similarities):.2f}%\n\n")
        
        # Distribution analysis
        self.stats_text.insert(tk.END, "Similarity Distribution:\n")
        
        ranges = [(90, 100), (80, 90), (70, 80), (60, 70), (50, 60), (0, 50)]
        for low, high in ranges:
            count = sum(1 for s in similarities if low <= s < high)
            percentage = (count / num_alignments * 100) if num_alignments > 0 else 0
            bar = "█" * int(percentage / 2)
            self.stats_text.insert(tk.END, f"  {low:3d}-{high:3d}%: {bar} ({count} regions, {percentage:.1f}%)\n")
        
        self.stats_text.insert(tk.END, "\n")
        
        # Top 10 alignments
        self.stats_text.insert(tk.END, "Top 10 Alignment Regions:\n")
        self.stats_text.insert(tk.END, f"{'Rank':<6} {'Score':<8} {'Similarity':<12} {'Region 1':<20} {'Region 2':<20}\n")
        self.stats_text.insert(tk.END, "-"*80 + "\n")
        
        for idx, result in enumerate(self.alignment_results[:10], 1):
            region1 = f"{result['window1_start']}-{result['window1_end']}"
            region2 = f"{result['window2_start']}-{result['window2_end']}"
            self.stats_text.insert(tk.END, 
                f"{idx:<6} {result['score']:<8} {result['similarity']:>6.2f}%     "
                f"{region1:<20} {region2:<20}\n")
    
    def export_results(self):
        """Export results to a text file"""
        if not self.alignment_results:
            messagebox.showwarning("Warning", "No results to export!")
            return
        
        filename = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
            initialfile=f"alignment_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
        )
        
        if filename:
            try:
                with open(filename, 'w') as f:
                    # Write header
                    f.write("="*100 + "\n")
                    f.write("GENOME LOCAL ALIGNMENT RESULTS\n")
                    f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                    f.write("="*100 + "\n\n")
                    
                    f.write(f"Genome 1: {self.genome1_name}\n")
                    f.write(f"Length: {len(self.genome1)} bp\n\n")
                    
                    f.write(f"Genome 2: {self.genome2_name}\n")
                    f.write(f"Length: {len(self.genome2)} bp\n\n")
                    
                    f.write(f"Parameters:\n")
                    f.write(f"  Match: {self.match_score}\n")
                    f.write(f"  Mismatch: {self.mismatch_score}\n")
                    f.write(f"  Gap: {self.gap_penalty}\n")
                    f.write(f"  Window Size: {self.window_size}\n")
                    f.write(f"  Step Size: {self.step_size}\n\n")
                    
                    # Write all results
                    for idx, result in enumerate(self.alignment_results, 1):
                        f.write(f"\nAlignment #{idx}\n")
                        f.write(f"{'-'*100}\n")
                        f.write(f"Region 1: {result['window1_start']}-{result['window1_end']} bp\n")
                        f.write(f"Region 2: {result['window2_start']}-{result['window2_end']} bp\n")
                        f.write(f"Score: {result['score']}\n")
                        f.write(f"Similarity: {result['similarity']:.2f}%\n")
                        f.write(f"Alignment length: {len(result['align1'])} positions\n\n")
                
                messagebox.showinfo("Success", f"Results exported to {filename}")
                self.status_var.set(f"Results exported to {filename}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to export: {str(e)}")
    
    def clear_results(self):
        """Clear all results"""
        self.results_text.delete(1.0, tk.END)
        self.viz_text.delete(1.0, tk.END)
        self.stats_text.delete(1.0, tk.END)
        self.alignment_results = []
        self.status_var.set("Results cleared")

def main():
    root = tk.Tk()
    app = GenomeLocalAlignment(root)
    root.mainloop()

if __name__ == "__main__":
    main()
