from pathlib import Path

# pyright: reportMissingModuleSource=false
import matplotlib.pyplot as plt
import pandas as pd


def main():
    results = pd.read_csv(
        "results/benchmark_results.csv"
    )

    Path("results").mkdir(exist_ok=True)

    plt.figure(figsize=(8, 5))

    plt.plot(
        results["length"],
        results["custom_time_seconds"],
        marker="o",
        label="Custom Python",
    )

    plt.plot(
        results["length"],
        results["biopython_time_seconds"],
        marker="o",
        label="BioPython",
    )

    plt.xlabel("Sequence length")
    plt.ylabel("Runtime in seconds")
    plt.title("Custom Aligner vs BioPython")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        "results/runtime_comparison.png",
        dpi=200,
    )

    plt.show()


if __name__ == "__main__":
    main()