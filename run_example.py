import sys
from pathlib import Path

try:
    import matplotlib.pyplot as plt  
except (ModuleNotFoundError, ImportError):  
    plt = None

try:
    from Bio.Align import substitution_matrices  
except (ModuleNotFoundError, ImportError):  
    substitution_matrices = None

sys.path.insert(0, str(Path(__file__).parent / "src"))

from alignment_tool import align_affine  


def ensure_results_dir():
    results_dir = Path("results")
    if results_dir.exists() and not results_dir.is_dir():
        results_dir.unlink()
    results_dir.mkdir(exist_ok=True)
    return results_dir


def plot_matrix(result, title, filename):
    if plt is None:
        raise ModuleNotFoundError("matplotlib is required to plot score matrices")
    plt.figure(figsize=(9, 6))
    plt.imshow(
        result.score_matrix,
        origin="lower",
        aspect="auto",
        cmap="viridis",
    )
    plt.colorbar(label="Alignment score")
    plt.xlabel("Position in sequence B")
    plt.ylabel("Position in sequence A")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(filename, dpi=200)
    plt.show()
    plt.close()


def main():
    results_dir = ensure_results_dir()

    if substitution_matrices is None:
        raise ModuleNotFoundError("biopython is required to load substitution matrices")

    matrix = substitution_matrices.load("BLOSUM62")

    sequence_a = "HEAGAWGHEE"
    sequence_b = "PAWHEAE"

    for mode in ["global", "local"]:
        result = align_affine(
            sequence_a,
            sequence_b,
            matrix,
            gap_open=10,
            gap_extend=0.5,
            mode=mode,
        )

        print("=" * 50)
        print(mode.upper())
        print("Score:", result.score)
        print(result.aligned_a)
        print(result.aligned_b)

        title = (
            "Needleman-Wunsch Score Matrix"
            if mode == "global"
            else "Smith-Waterman Score Matrix"
        )

        filename = (
            results_dir / f"{mode}_score_matrix.png"
        )

        plot_matrix(result, title, filename)


if __name__ == "__main__":
    main()