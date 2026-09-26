# Sequence Alignment Tool

This project implements a custom protein sequence-alignment tool in Python using dynamic programming. It supports both the Needleman–Wunsch algorithm for global alignment and the Smith–Waterman algorithm for local alignment. The implementation is written from scratch rather than relying on BioPython for the main alignment calculations.

The tool uses the BLOSUM62 substitution matrix to score amino-acid matches and substitutions. It also supports affine gap penalties, which provide separate costs for opening and extending a gap. Three dynamic-programming matrices are maintained: one for aligned residues and two for gaps in either sequence. This approach models biological insertions and deletions more realistically than a constant gap penalty.

For each alignment, the program calculates the scoring matrices, performs traceback, and returns the aligned sequences and final score. Needleman–Wunsch traceback begins at the bottom-right cell because it aligns both complete sequences. Smith–Waterman traceback begins at the highest-scoring cell and stops when the score reaches zero, identifying the best local matching region.

The project includes visualizations of the dynamic-programming matrices using matplotlib. These plots help demonstrate how alignment scores develop across the two sequences. The custom implementation is also tested against BioPython’s optimized PairwiseAligner to verify that the calculated alignment scores are correct for both global and local alignment modes.

Runtime benchmarking is included to compare the pure Python implementation with BioPython. BioPython is expected to be faster because its alignment routines are optimized and implemented using compiled code. Both implementations have (O(nm)) time complexity, where (n) and (m) are the sequence lengths, while storing all affine-gap matrices requires (O(nm)) memory.




### Needleman-Wunsch and Smith-Waterman with affine gap penalties

The alignment score is computed using dynamic programming with three matrices: $M$, $X$, and $Y$.

- $M(i,j)$ represents a match or mismatch between residues $a_i$ and $b_j$
- $X(i,j)$ represents a gap in sequence $B$
- $Y(i,j)$ represents a gap in sequence $A$

For affine gap penalties, the recurrence is:

$$
M(i,j) = \max\{M(i-1,j-1), X(i-1,j-1), Y(i-1,j-1)\} + s(a_i,b_j)
$$

$$
X(i,j) = \max\{M(i-1,j)-g_{open}, X(i-1,j)-g_{extend}\}
$$

$$
Y(i,j) = \max\{M(i,j-1)-g_{open}, Y(i,j-1)-g_{extend}\}
$$

where $s(a_i,b_j)$ is the substitution score and $g_{open}$ and $g_{extend}$ are the gap-opening and gap-extension penalties.

For global alignment, the optimal score is selected from the final cell:

$$
\text{score}_{global} = \max\{M(n,m), X(n,m), Y(n,m)\}
$$

For local alignment, the recurrence is modified with a reset to zero:

$$
M(i,j) = \max\{0,\; M(i-1,j-1),\; X(i-1,j-1),\; Y(i-1,j-1)\} + s(a_i,b_j)
$$

$$
X(i,j) = \max\{0,\; M(i-1,j)-g_{open},\; X(i-1,j)-g_{extend}\}
$$

$$
Y(i,j) = \max\{0,\; M(i,j-1)-g_{open},\; Y(i,j-1)-g_{extend}\}
$$

and the local alignment score is:

$$
\text{score}_{local} = \max_{i,j}\{M(i,j), X(i,j), Y(i,j)\}
$$

This formulation allows the algorithm to handle both global sequence alignment (Needleman-Wunsch) and local sequence alignment (Smith-Waterman) while correctly modeling affine gap penalties.




# global / local affine gap recurrence
M[i, j] = max(M[i-1, j-1], X[i-1, j-1], Y[i-1, j-1]) + s[a[i-1], b[j-1]]
X[i, j] = max(M[i-1, j] - gap_open, X[i-1, j] - gap_extend)
Y[i, j] = max(M[i, j-1] - gap_open, Y[i, j-1] - gap_extend)

# local alignment reset
M[i, j] = max(0.0, M[i-1, j-1], X[i-1, j-1], Y[i-1, j-1]) + s[a[i-1], b[j-1]]
X[i, j] = max(0.0, M[i-1, j] - gap_open, X[i-1, j] - gap_extend)
Y[i, j] = max(0.0, M[i, j-1] - gap_open, Y[i, j-1] - gap_extend)
``

## Installation

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt