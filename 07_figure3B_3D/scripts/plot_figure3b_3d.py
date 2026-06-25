from __future__ import annotations

import argparse
import csv
import random
import re
import statistics
from pathlib import Path

import matplotlib.pyplot as plt

from figure_metrics import model_metrics, read_measured_activity, read_rates

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
DEFAULT_MMC4 = DATA / "mmc_files" / "1-s2.0-S2589004223028249-mmc4.xlsx"

EXPECTED = {1: 7, 2: 21, 3: 35, 4: 35, 5: 21, 6: 7, 7: 1}
FIGURE3D_FOLDERS = {
    1: "1TF_full",
    2: "2TF_medium",
    3: "3TF_medium",
    4: "4TF_medium",
    5: "5TF_medium",
    6: "6TF_full",
    7: "7TF_full",
}


def quartiles(values: list[float]) -> tuple[float, float]:
    values = sorted(values)
    if len(values) < 4:
        return min(values), max(values)
    q1, _, q3 = statistics.quantiles(values, n=4, method="inclusive")
    return q1, q3


def combination_name(path: Path, tf_count: int) -> str:
    stem = path.name.removesuffix("_rate_invert.txt")
    if "_selfoff_" in stem:
        if tf_count == 1:
            match = re.search(r"_([1-7])_scramble$", stem)
            index_to_tf = {
                "1": "ATF1",
                "2": "ATF4",
                "3": "ATF7",
                "4": "CREB1",
                "5": "CREB3",
                "6": "CREB5",
                "7": "CREM",
            }
            return index_to_tf.get(match.group(1) if match else "", stem)
        match = re.search(r"_comb_\d+_(.+)$", stem)
        return match.group(1).replace("_", "-") if match else stem
    if tf_count == 1:
        match = re.search(r"1TF_(.+?)_(?:rep\d+_)?seed", stem)
    elif tf_count == 7:
        return "CREB1-CREB3-CREB5-CREM-ATF1-ATF4-ATF7"
    else:
        match = re.search(rf"{tf_count}TF_(.+?)_(?:medium_)?(?:full_)?seed", stem)
    return match.group(1) if match else stem


def files_for_figure(root: Path, figure: str, tf_count: int) -> list[Path]:
    if figure == "3B":
        folder = root / "figure3b_selfcompetition_off" / f"{tf_count}TF"
        if tf_count == 1:
            # Exclude the extra named CREB1 file; Figure 3B uses the seven numbered 1TF models.
            return sorted(folder.glob("CRE_1TF_selfoff_singlehit_fullanneal_[1-7]_scramble_rate_invert.txt"))
        return sorted(folder.glob("*_comb_*_rate_invert.txt"))
    folder = root / "figure3d_selfcompetition_on" / FIGURE3D_FOLDERS[tf_count]
    return sorted(folder.glob("*_rate_invert.txt"))


def calculate(mmc4: Path, figure: str) -> list[dict[str, object]]:
    measured = read_measured_activity(mmc4)
    data_root = DATA
    results: list[dict[str, object]] = []
    for tf_count in range(1, 8):
        files = files_for_figure(data_root, figure, tf_count)
        if len(files) != EXPECTED[tf_count]:
            raise ValueError(f"{figure} {tf_count}TF: expected {EXPECTED[tf_count]}, found {len(files)}")
        for path in files:
            r_value, rmse, sample_count = model_metrics(read_rates(path), measured)
            results.append(
                {
                    "figure": figure,
                    "self_competition": "off" if figure == "3B" else "on",
                    "TF_count": tf_count,
                    "combination": combination_name(path, tf_count),
                    "R": r_value,
                    "R_squared": r_value**2,
                    "RMSE": rmse,
                    "N": sample_count,
                    "file": str(path.relative_to(ROOT)),
                }
            )
    return results


def save_metrics(results: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)


def draw(results: list[dict[str, object]], figure: str, output: Path, ymin: float, ymax: float, dpi: int, show: bool) -> None:
    grouped = {
        count: [float(row["R"]) for row in results if int(row["TF_count"]) == count]
        for count in range(1, 8)
    }
    title = (
        "Figure 3B | Single-hit fitting without self-competition"
        if figure == "3B"
        else "Figure 3D | Single-hit fitting with self-competition"
    )
    color = "#4D86A8" if figure == "3B" else "#67A695"

    plt.rcParams.update({"font.family": "Arial", "font.size": 10, "savefig.facecolor": "white"})
    fig, ax = plt.subplots(figsize=(8.2, 5.8))
    boxes = ax.boxplot(
        [grouped[count] for count in range(1, 8)],
        positions=range(1, 8),
        widths=0.58,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#202428", "linewidth": 1.5},
        whiskerprops={"color": "#50575C", "linewidth": 1.0},
        capprops={"color": "#50575C", "linewidth": 1.0},
        boxprops={"edgecolor": "#356D62" if figure == "3D" else "#3A657D", "linewidth": 1.0},
    )
    for box in boxes["boxes"]:
        box.set_facecolor(color)
        box.set_alpha(0.76)

    rng = random.Random(2024)
    for count, values in grouped.items():
        q1, q3 = quartiles(values)
        iqr = q3 - q1
        outliers = [v for v in values if v < q1 - 1.5 * iqr or v > q3 + 1.5 * iqr]
        if outliers:
            ax.scatter(
                [count + rng.uniform(-0.05, 0.05) for _ in outliers],
                outliers,
                s=28,
                color="#D4573F",
                edgecolor="white",
                linewidth=0.45,
                zorder=4,
            )
        ax.scatter(
            count,
            statistics.fmean(values),
            s=38,
            color="black",
            edgecolor="white",
            linewidth=0.55,
            zorder=5,
            label="Group mean" if count == 1 else None,
        )

    ax.axhline(0.8, color="#B14D47", linestyle=(0, (4, 3)), label="R = 0.80")
    ax.axhline(0, color="#777777", linewidth=0.7)
    ax.set_xlim(0.5, 7.5)
    ax.set_ylim(ymin, ymax)
    ax.set_xticks(range(1, 8))
    ax.set_xlabel("Number of ATF/CREB family TFs (n)")
    ax.set_ylabel("Pearson's R")
    ax.set_title(title, fontsize=13, fontweight="bold", pad=14)
    ax.grid(axis="y", color="#D8DDE0", linewidth=0.65, alpha=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, loc="lower right")
    for count, values in grouped.items():
        ax.text(count, ymin + 0.025, f"n={len(values)}", ha="center", fontsize=8)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)


def run_one(figure: str, mmc4: Path, dpi: int, show: bool) -> None:
    output = ROOT / "outputs" / (
        "Figure3B_selfcompetition_off.png" if figure == "3B" else "Figure3D_selfcompetition_on.png"
    )
    ymin, ymax = (-0.2, 1.0)
    results = calculate(mmc4, figure)
    metrics = output.with_name(f"{output.stem}_metrics.csv")
    save_metrics(results, metrics)
    draw(results, figure, output, ymin, ymax, dpi, show)
    print(f"{figure}: saved {output}")
    for count in range(1, 8):
        values = [float(row["R"]) for row in results if int(row["TF_count"]) == count]
        print(f"{figure} {count}TF: n={len(values):2d}, mean={statistics.fmean(values):.4f}, median={statistics.median(values):.4f}, max={max(values):.4f}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Reproduce Figure 3B and Figure 3D.")
    parser.add_argument("--figure", choices=["3B", "3D", "both"], default="both")
    parser.add_argument("--mmc4", type=Path, default=DEFAULT_MMC4)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    figures = ["3B", "3D"] if args.figure == "both" else [args.figure]
    for figure in figures:
        run_one(figure, args.mmc4, args.dpi, args.show)


if __name__ == "__main__":
    main()
