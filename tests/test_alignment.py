import sys
from pathlib import Path

from Bio.Align import PairwiseAligner  
from Bio.Align import substitution_matrices  

sys.path.insert(
    0,
    str(Path(__file__).parents[1] / "src"),
)

import importlib
    

try:
    alignment_mod = importlib.import_module("alignment_tool")
except Exception:
    
    sys.path.insert(
        0,
        str(Path(__file__).parents[1] / "src"),
    )
    alignment_mod = importlib.import_module("alignment_tool")

align_affine = alignment_mod.align_affine


def make_reference_aligner(mode, matrix):
    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    return aligner


def test_global_score_matches_biopython():
    matrix = substitution_matrices.load("BLOSUM62")
    sequence_a = "HEAGAWGHEE"
    sequence_b = "PAWHEAE"

    custom = align_affine(
        sequence_a,
        sequence_b,
        matrix,
        10,
        0.5,
        "global",
    )

    reference = make_reference_aligner("global", matrix)
    expected = reference.score(sequence_a, sequence_b)

    assert abs(custom.score - expected) < 1e-6


def test_local_score_matches_biopython():
    matrix = substitution_matrices.load("BLOSUM62")
    sequence_a = "HEAGAWGHEE"
    sequence_b = "PAWHEAE"

    custom = align_affine(
        sequence_a,
        sequence_b,
        matrix,
        10,
        0.5,
        "local",
    )

    reference = make_reference_aligner("local", matrix)
    expected = reference.score(sequence_a, sequence_b)

    assert abs(custom.score - expected) < 1e-6


def test_alignment_lengths_are_equal():
    matrix = substitution_matrices.load("BLOSUM62")

    result = align_affine(
        "HEAGAWGHEE",
        "PAWHEAE",
        matrix,
        10,
        0.5,
        "local",
    )

    assert len(result.aligned_a) == len(result.aligned_b)