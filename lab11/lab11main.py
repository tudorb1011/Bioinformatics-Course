import sys
import numpy as np
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout,
                             QHBoxLayout, QLabel, QLineEdit, QPushButton,
                             QTextEdit, QGroupBox, QCheckBox, QGridLayout)
from PyQt5.QtCore import Qt
import matplotlib

matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure


class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)


class DNAAlignmentApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("DNA Sequence Alignment - Needleman-Wunsch Algorithm")
        self.setGeometry(100, 100, 1400, 800)

        # Default sequences
        self.default_s1 = "ACCGTGAAGCCAATAC"
        self.default_s2 = "AGCGTGCAGCCAATAC"

        self.init_ui()

    def init_ui(self):
        # Main widget and layout
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QHBoxLayout()
        main_widget.setLayout(main_layout)

        # Left panel
        left_panel = QVBoxLayout()
        left_panel.addWidget(self.create_sequence_group())
        left_panel.addWidget(self.create_parameters_group())
        left_panel.addStretch()

        # Middle panel - Visualizations
        middle_panel = QVBoxLayout()
        self.matrix_canvas = MplCanvas(self, width=6, height=5, dpi=100)
        matrix_group = QGroupBox("Graphic representation of the alignment matrix")
        matrix_layout = QVBoxLayout()
        matrix_layout.addWidget(self.matrix_canvas)
        matrix_group.setLayout(matrix_layout)
        middle_panel.addWidget(matrix_group)

        self.traceback_canvas = MplCanvas(self, width=6, height=3, dpi=100)
        traceback_group = QGroupBox("Traceback path deviation from optimal alignment (diagonal)")
        traceback_layout = QVBoxLayout()
        traceback_layout.addWidget(self.traceback_canvas)
        traceback_group.setLayout(traceback_layout)
        middle_panel.addWidget(traceback_group)

        # Right panel - Results
        right_panel = QVBoxLayout()
        results_group = QGroupBox("Show Alignment:")
        results_layout = QVBoxLayout()
        self.results_text = QTextEdit()
        self.results_text.setReadOnly(True)
        self.results_text.setFontFamily("Courier")
        self.results_text.setMinimumWidth(400)
        results_layout.addWidget(self.results_text)
        results_group.setLayout(results_layout)
        right_panel.addWidget(results_group)

        # Add all panels to main layout
        left_container = QWidget()
        left_container.setLayout(left_panel)
        left_container.setMaximumWidth(350)

        middle_container = QWidget()
        middle_container.setLayout(middle_panel)

        right_container = QWidget()
        right_container.setLayout(right_panel)

        main_layout.addWidget(left_container)
        main_layout.addWidget(middle_container)
        main_layout.addWidget(right_container)

    def create_sequence_group(self):
        group = QGroupBox("Sequences")
        layout = QGridLayout()

        layout.addWidget(QLabel("Sq 1 ="), 0, 0)
        self.seq1_entry = QLineEdit(self.default_s1)
        layout.addWidget(self.seq1_entry, 0, 1)

        layout.addWidget(QLabel("Sq 2 ="), 1, 0)
        self.seq2_entry = QLineEdit(self.default_s2)
        layout.addWidget(self.seq2_entry, 1, 1)

        group.setLayout(layout)
        return group

    def create_parameters_group(self):
        group = QGroupBox("Parameters")
        layout = QGridLayout()

        layout.addWidget(QLabel("Gap ="), 0, 0)
        self.gap_entry = QLineEdit("0")
        self.gap_entry.setMaximumWidth(80)
        layout.addWidget(self.gap_entry, 0, 1)

        layout.addWidget(QLabel("Match ="), 1, 0)
        self.match_entry = QLineEdit("1")
        self.match_entry.setMaximumWidth(80)
        layout.addWidget(self.match_entry, 1, 1)

        layout.addWidget(QLabel("MMach ="), 2, 0)
        self.mismatch_entry = QLineEdit("-1")
        self.mismatch_entry.setMaximumWidth(80)
        layout.addWidget(self.mismatch_entry, 2, 1)

        # Options
        self.plot_grid_check = QCheckBox("Plot grid")
        self.plot_grid_check.setChecked(True)
        layout.addWidget(self.plot_grid_check, 3, 0, 1, 2)

        # Align button
        align_btn = QPushButton("Align")
        align_btn.setMinimumHeight(50)
        align_btn.clicked.connect(self.perform_alignment)
        layout.addWidget(align_btn, 4, 0, 1, 2)

        group.setLayout(layout)
        return group

    def needleman_wunsch(self, seq1, seq2, match_score, mismatch_score, gap_penalty):
        """Implements the Needleman-Wunsch algorithm for global sequence alignment"""
        n, m = len(seq1), len(seq2)

        # Initialize scoring matrix
        score_matrix = np.zeros((n + 1, m + 1))

        # Initialize first row and column
        for i in range(n + 1):
            score_matrix[i][0] = gap_penalty * i
        for j in range(m + 1):
            score_matrix[0][j] = gap_penalty * j

        # Fill the scoring matrix
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                match = score_matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)
                delete = score_matrix[i - 1][j] + gap_penalty
                insert = score_matrix[i][j - 1] + gap_penalty
                score_matrix[i][j] = max(match, delete, insert)

        # Traceback
        align1, align2 = "", ""
        i, j = n, m
        traceback_path = [(i, j)]

        while i > 0 or j > 0:
            current_score = score_matrix[i][j]

            if i > 0 and j > 0:
                diagonal_score = score_matrix[i - 1][j - 1]
                match_mismatch = match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score

                if current_score == diagonal_score + match_mismatch:
                    align1 = seq1[i - 1] + align1
                    align2 = seq2[j - 1] + align2
                    i -= 1
                    j -= 1
                    traceback_path.append((i, j))
                    continue

            if i > 0 and current_score == score_matrix[i - 1][j] + gap_penalty:
                align1 = seq1[i - 1] + align1
                align2 = "-" + align2
                i -= 1
                traceback_path.append((i, j))
            elif j > 0 and current_score == score_matrix[i][j - 1] + gap_penalty:
                align1 = "-" + align1
                align2 = seq2[j - 1] + align2
                j -= 1
                traceback_path.append((i, j))
            else:
                break

        traceback_path.reverse()
        return score_matrix, align1, align2, traceback_path

    def calculate_statistics(self, align1, align2):
        """Calculate alignment statistics"""
        matches = sum(1 for a, b in zip(align1, align2) if a == b and a != '-')
        length = len(align1)
        similarity = (matches / length * 100) if length > 0 else 0
        return matches, length, similarity

    def visualize_matrix(self, score_matrix, traceback_path):
        """Visualize the scoring matrix with heatmap"""
        self.matrix_canvas.axes.clear()

        # Create heatmap
        im = self.matrix_canvas.axes.imshow(score_matrix, cmap='RdPu', aspect='auto', interpolation='nearest')

        # Plot traceback path
        if traceback_path:
            path_i = [p[0] for p in traceback_path]
            path_j = [p[1] for p in traceback_path]
            self.matrix_canvas.axes.plot(path_j, path_i, 'b-', linewidth=2, alpha=0.7)

        if self.plot_grid_check.isChecked():
            self.matrix_canvas.axes.set_xticks(np.arange(score_matrix.shape[1]))
            self.matrix_canvas.axes.set_yticks(np.arange(score_matrix.shape[0]))
            self.matrix_canvas.axes.grid(True, which='both', color='white', linewidth=0.5)

        self.matrix_canvas.figure.tight_layout()
        self.matrix_canvas.draw()

    def visualize_traceback(self, traceback_path, seq1_len, seq2_len):
        """Visualize traceback path deviation from diagonal"""
        self.traceback_canvas.axes.clear()

        # Create grid
        grid = np.ones((seq1_len + 1, seq2_len + 1)) * 0.9

        # Mark traceback path
        for i, j in traceback_path:
            grid[i][j] = 0.2

        self.traceback_canvas.axes.imshow(grid, cmap='RdYlGn_r', aspect='auto', interpolation='nearest')

        if self.plot_grid_check.isChecked():
            self.traceback_canvas.axes.set_xticks(np.arange(seq2_len + 1))
            self.traceback_canvas.axes.set_yticks(np.arange(seq1_len + 1))
            self.traceback_canvas.axes.grid(True, which='both', color='white', linewidth=0.5)

        self.traceback_canvas.figure.tight_layout()
        self.traceback_canvas.draw()

    def format_alignment_display(self, align1, align2):
        """Format alignment with match indicators"""
        match_line = ""
        for a, b in zip(align1, align2):
            if a == b and a != '-':
                match_line += "|"
            else:
                match_line += " "
        return match_line

    def perform_alignment(self):
        """Perform the alignment and display results"""
        try:
            # Get input values
            seq1 = self.seq1_entry.text().upper().strip()
            seq2 = self.seq2_entry.text().upper().strip()
            gap_penalty = int(self.gap_entry.text())
            match_score = int(self.match_entry.text())
            mismatch_score = int(self.mismatch_entry.text())

            # Validate sequences
            if not seq1 or not seq2:
                self.results_text.setText("Error: Please enter both sequences")
                return

            # Perform alignment
            score_matrix, align1, align2, traceback_path = self.needleman_wunsch(
                seq1, seq2, match_score, mismatch_score, gap_penalty
            )

            # Calculate statistics
            matches, length, similarity = self.calculate_statistics(align1, align2)

            # Display results
            match_line = self.format_alignment_display(align1, align2)

            results = []
            results.append(f"A-{align1}")
            results.append(f"  {match_line}")
            results.append(f"AG-{align2}")
            results.append("")
            results.append(f"Maches = {matches}")
            results.append(f"Length = {length}")
            results.append("")
            results.append(f"Similarity = {similarity:.0f} %")
            results.append("")
            results.append(f"Tracing back: M[{len(seq1)},{len(seq2)}]")
            results.append("")
            results.append("Scoring Matrix:")
            results.append("    " + "  ".join([" "] + list(seq2)))

            for i in range(min(len(seq1) + 1, 10)):
                row_label = " " if i == 0 else seq1[i - 1]
                row_values = "  ".join([f"{int(score_matrix[i][j]):2d}" for j in range(min(len(seq2) + 1, 15))])
                results.append(f" {row_label}  {row_values}")

            self.results_text.setText("\n".join(results))

            # Visualize
            self.visualize_matrix(score_matrix, traceback_path)
            self.visualize_traceback(traceback_path, len(seq1), len(seq2))

        except Exception as e:
            self.results_text.setText(f"Error: {str(e)}")
            import traceback
            self.results_text.append("\n" + traceback.format_exc())


def main():
    app = QApplication(sys.argv)
    window = DNAAlignmentApp()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()