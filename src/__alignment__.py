from dataclasses import dataclass
from typing import Dict

import numpy as np  # pyright: ignore[reportMissingImports]


@dataclass
class AlignmentResult:
    aligned_a: str
    aligned_b: str
    score: float
    score_matrix: np.ndarray
    matrices: Dict[str, np.ndarray]


def _best(values):
    index = max(range(len(values)), key=lambda position: values[position])
    return values[index], index


def align_affine(
    sequence_a: str,
    sequence_b: str,
    substitution_matrix,
    gap_open: float,
    gap_extend: float,
    mode: str = "global",
) -> AlignmentResult:
    """
    Align two protein sequences using affine gap penalties.

    mode:
        "global" -> Needleman-Wunsch
        "local"  -> Smith-Waterman

    A gap of length k costs:
        gap_open + (k - 1) * gap_extend
    """

    if mode not in {"global", "local"}:
        raise ValueError("mode must be 'global' or 'local'")

    if gap_open < 0 or gap_extend < 0:
        raise ValueError("Gap penalties must be non-negative")

    sequence_a = sequence_a.upper()
    sequence_b = sequence_b.upper()

    n = len(sequence_a)
    m = len(sequence_b)
    negative_infinity = -np.inf

    # M: residue aligned with residue
    # X: residue from A aligned with a gap
    # Y: gap aligned with residue from B
    matrix_m = np.full((n + 1, m + 1), negative_infinity)
    matrix_x = np.full((n + 1, m + 1), negative_infinity)
    matrix_y = np.full((n + 1, m + 1), negative_infinity)

    pointer_m = np.full((n + 1, m + 1), "", dtype=object)
    pointer_x = np.full((n + 1, m + 1), "", dtype=object)
    pointer_y = np.full((n + 1, m + 1), "", dtype=object)

    matrix_m[0, 0] = 0.0

    if mode == "local":
        matrix_m[:, :] = 0.0
        matrix_x[:, :] = 0.0
        matrix_y[:, :] = 0.0
    else:
        for i in range(1, n + 1):
            matrix_x[i, 0] = -gap_open - (i - 1) * gap_extend
            pointer_x[i, 0] = "X" if i > 1 else ""

        for j in range(1, m + 1):
            matrix_y[0, j] = -gap_open - (j - 1) * gap_extend
            pointer_y[0, j] = "Y" if j > 1 else ""

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            substitution_score = float(
                substitution_matrix[
                    sequence_a[i - 1],
                    sequence_b[j - 1],
                ]
            )

            diagonal_values = [
                matrix_m[i - 1, j - 1],
                matrix_x[i - 1, j - 1],
                matrix_y[i - 1, j - 1],
            ]

            previous_score, previous_index = _best(diagonal_values)
            candidate_m = previous_score + substitution_score

            if mode == "local" and candidate_m < 0:
                matrix_m[i, j] = 0.0
                pointer_m[i, j] = ""
            else:
                matrix_m[i, j] = candidate_m
                pointer_m[i, j] = ("M", "X", "Y")[previous_index]

            opening_x = matrix_m[i - 1, j] - gap_open
            extending_x = matrix_x[i - 1, j] - gap_extend

            if mode == "local":
                value, index = _best([0.0, opening_x, extending_x])
                matrix_x[i, j] = value
                pointer_x[i, j] = ("", "M", "X")[index]
            elif opening_x >= extending_x:
                matrix_x[i, j] = opening_x
                pointer_x[i, j] = "M"
            else:
                matrix_x[i, j] = extending_x
                pointer_x[i, j] = "X"

            opening_y = matrix_m[i, j - 1] - gap_open
            extending_y = matrix_y[i, j - 1] - gap_extend

            if mode == "local":
                value, index = _best([0.0, opening_y, extending_y])
                matrix_y[i, j] = value
                pointer_y[i, j] = ("", "M", "Y")[index]
            elif opening_y >= extending_y:
                matrix_y[i, j] = opening_y
                pointer_y[i, j] = "M"
            else:
                matrix_y[i, j] = extending_y
                pointer_y[i, j] = "Y"

    if mode == "global":
        final_values = [matrix_m[n, m], matrix_x[n, m], matrix_y[n, m]]
        score, state_index = _best(final_values)
        state = ("M", "X", "Y")[state_index]
        i, j = n, m
    else:
        score_matrix = np.maximum.reduce([matrix_m, matrix_x, matrix_y])
        i, j = np.unravel_index(np.argmax(score_matrix), score_matrix.shape)
        score = score_matrix[i, j]

        final_values = [matrix_m[i, j], matrix_x[i, j], matrix_y[i, j]]
        _, state_index = _best(final_values)
        state = ("M", "X", "Y")[state_index]

    aligned_a = []
    aligned_b = []

    while i > 0 or j > 0:
        if mode == "local":
            current_score = {
                "M": matrix_m[i, j],
                "X": matrix_x[i, j],
                "Y": matrix_y[i, j],
            }[state]

            if current_score <= 0:
                break

        if state == "M":
            if i == 0 or j == 0:
                break

            aligned_a.append(sequence_a[i - 1])
            aligned_b.append(sequence_b[j - 1])
            state = pointer_m[i, j]
            i -= 1
            j -= 1

        elif state == "X":
            if i == 0:
                break

            aligned_a.append(sequence_a[i - 1])
            aligned_b.append("-")
            previous_state = pointer_x[i, j]
            i -= 1

            if previous_state == "":
                break

            state = previous_state

        elif state == "Y":
            if j == 0:
                break

            aligned_a.append("-")
            aligned_b.append(sequence_b[j - 1])
            previous_state = pointer_y[i, j]
            j -= 1

            if previous_state == "":
                break

            state = previous_state

        else:
            raise RuntimeError(f"Invalid traceback state: {state}")

    aligned_a.reverse()
    aligned_b.reverse()

    score_matrix = np.maximum.reduce([matrix_m, matrix_x, matrix_y])

    return AlignmentResult(
        aligned_a="".join(aligned_a),
        aligned_b="".join(aligned_b),
        score=float(score),
        score_matrix=score_matrix,
        matrices={
            "M": matrix_m,
            "X": matrix_x,
            "Y": matrix_y,
        },
    )