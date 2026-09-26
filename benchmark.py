import random
import csv
import sys
import time
from pathlib import Path

from Bio.Align import PairwiseAligner  
from Bio.Align import substitution_matrices  

PROJECT_ROOT = Path(__file__).resolve().parent
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from alignment_tool import align_affine  


def ensure_results_dir():
    results_dir = Path("results")
    if results_dir.exists() and not results_dir.is_dir():
        results_dir.unlink()
    results_dir.mkdir(exist_ok=True)
    return results_dir


def random_protein(length):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    return "".join(random.choice(amino_acids) for _ in range(length))


def main():
    matrix = substitution_matrices.load("BLOSUM62")

    reference = PairwiseAligner()
    reference.mode = "global"
    reference.substitution_matrix = matrix
    reference.open_gap_score = -10
    reference.extend_gap_score = -0.5

    rows = []

    for length in [50, 100, 200, 300, 500]:
        sequence_a = random_protein(length)
        sequence_b = random_protein(length)

        start = time.perf_counter()

        custom = align_affine(
            sequence_a,
            sequence_b,
            matrix,
            10,
            0.5,
            "global",
        )

        custom_time = time.perf_counter() - start

        start = time.perf_counter()
        reference_score = reference.score(sequence_a, sequence_b)
        biopython_time = time.perf_counter() - start

        rows.append({
            "length": length,
            "custom_score": custom.score,
            "biopython_score": reference_score,
            "scores_match": abs(
                custom.score - reference_score
            ) < 1e-6,
            "custom_time_seconds": custom_time,
            "biopython_time_seconds": biopython_time,
        })

    results_dir = ensure_results_dir()
    output_path = results_dir / "benchmark_results.csv"

    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    headers = list(rows[0].keys())
    print(" ".join(f"{header:>22}" for header in headers))
    for row in rows:
        print(
            " ".join(
                f"{(f'{value:.6f}' if isinstance(value, float) else value):>22}"
                for value in row.values()
            )
        )


if __name__ == "__main__":
    main()