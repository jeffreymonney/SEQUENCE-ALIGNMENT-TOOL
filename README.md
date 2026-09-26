# Sequence Alignment Tool

This project implements a custom protein sequence-alignment tool in Python using dynamic programming. It supports both the Needleman–Wunsch algorithm for global alignment and the Smith–Waterman algorithm for local alignment. The implementation is written from scratch rather than relying on BioPython for the main alignment calculations.

The tool uses the BLOSUM62 substitution matrix to score amino-acid matches and substitutions. It also supports affine gap penalties, which provide separate costs for opening and extending a gap. Three dynamic-programming matrices are maintained: one for aligned residues and two for gaps in either sequence. This approach models biological insertions and deletions more realistically than a constant gap penalty.

For each alignment, the program calculates the scoring matrices, performs traceback, and returns the aligned sequences and final score. Needleman–Wunsch traceback begins at the bottom-right cell because it aligns both complete sequences. Smith–Waterman traceback begins at the highest-scoring cell and stops when the score reaches zero, identifying the best local matching region.

The project includes visualizations of the dynamic-programming matrices using matplotlib. These plots help demonstrate how alignment scores develop across the two sequences. The custom implementation is also tested against BioPython’s optimized PairwiseAligner to verify that the calculated alignment scores are correct for both global and local alignment modes.

Runtime benchmarking is included to compare the pure Python implementation with BioPython. BioPython is expected to be faster because its alignment routines are optimized and implemented using compiled code. Both implementations have (O(nm)) time complexity, where (n) and (m) are the sequence lengths, while storing all affine-gap matrices requires (O(nm)) memory.

## Installation

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt