from __future__ import annotations

import argparse
import csv
import re
import statistics
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from figure_metrics import model_metrics, read_measured_activity, read_rates

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
DEFAULT_MMC4 = DATA / "mmc_files" / "1-s2.0-S2589004223028249-mmc4.xlsx"
DEFAULT_OUTPUT = ROOT / "outputs" / "Figure3A_TF_addition.png"

BASELINE_DIR = DATA / "rate_invert" / "creb1_self"
NON_FAMILY_DIR = DATA / "rate_invert" / "non_family"
FAMILY_DIR = DATA / "rate_invert" / "family"
EXCLUDED_RUNS = DATA / "Figure3A_excluded_nan_runs.txt"

NON_FAMILY_ORDER = [
    "DLX2", "HOXA10", "HOXA11", "HOXA5", "HOXA6", "HOXA9", "HOXB9",
    "HOXD13", "NKX22", "POU3F2", "SOX4", "TCF3", "YY1", "ZIC2",
]
FAMILY_ORDER = ["ATF1", "ATF4", "ATF7", "CREB3", "CREB5", "CREM"]
DISPLAY_NAMES = {"NKX22": "NKX2-2"}
COLORS = {"baseline": "#666D73", "non-family": "#D49A45", "family": "#4B86A8"}


def files_for(directory: Path, tf_name: str | None = None) -> list[Path]:
    if tf_name is None:
        return sorted(directory.glob("*CREB1_scramble[1-5]_rate_invert.txt"))
    return sorted(directory.glob(f"*_CREB1_{tf_name}_scramble*_rate_invert.txt"))


def calculate(mmc4: Path) -> list[dict[str, object]]:
    measured = read_measured_activity(mmc4)
    nan_runs = set()
    if EXCLUDED_RUNS.exists():
        for line in EXCLUDED_RUNS.read_text(encoding="utf-8").splitlines():
            fields = line.split("\t", maxsplit=1)
            if len(fields) == 2 and fields[0] == "NAN":
                nan_runs.add(Path(fields[1]).stem)

    groups = [("CREB1 self", "baseline", files_for(BASELINE_DIR))]
    groups += [(tf, "non-family", files_for(NON_FAMILY_DIR, tf)) for tf in NON_FAMILY_ORDER]
    groups += [(tf, "family", files_for(FAMILY_DIR, tf)) for tf in FAMILY_ORDER]

    results: list[dict[str, object]] = []
    for tf_name, group, files in groups:
        if not files:
            raise ValueError(f"No rate/invert files found for {tf_name}")
        for path in files:
            r_value, rmse, sample_count = model_metrics(read_rates(path), measured)
            match = re.search(r"scramble_?(\d+)", path.name)
            run_stem = path.name.removesuffix("_rate_invert.txt")
            results.append(
                {
                    "group": group,
                    "TF": tf_name,
                    "repeat": int(match.group(1)) if match else "",
                    "annealing_status": "s=-nan" if run_stem in nan_runs else "valid",
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


def draw(results: list[dict[str, object]], output: Path, ymin: float, ymax: float, dpi: int, show: bool) -> None:
    names = ["CREB1 self", *NON_FAMILY_ORDER, *FAMILY_ORDER]
    group_for = {"CREB1 self": "baseline", **{n: "non-family" for n in NON_FAMILY_ORDER}, **{n: "family" for n in FAMILY_ORDER}}
    grouped = {
        name: [float(row["R"]) for row in results if row["TF"] == name and row["annealing_status"] == "valid"]
        for name in names
    }
    baseline_best = max(grouped["CREB1 self"])
    positions = list(range(1, len(names) + 1))

    plt.rcParams.update({"font.family": "Arial", "font.size": 9, "savefig.facecolor": "white"})
    fig, ax = plt.subplots(figsize=(15.5, 6.8))
    boxes = ax.boxplot(
        [grouped[name] for name in names],
        positions=positions,
        widths=0.62,
        patch_artist=True,
        showfliers=True,
        medianprops={"color": "#1F2529", "linewidth": 1.35},
        whiskerprops={"color": "#50575C", "linewidth": 0.9},
        capprops={"color": "#50575C", "linewidth": 0.9},
        flierprops={"marker": "o", "markerfacecolor": "#202428", "markeredgecolor": "#202428", "markersize": 3.8},
    )
    for box, name in zip(boxes["boxes"], names):
        box.set_facecolor(COLORS[group_for[name]])
        box.set_edgecolor("#485157")
        box.set_alpha(0.78)

    for pos, name in zip(positions, names):
        nan_values = [float(row["R"]) for row in results if row["TF"] == name and row["annealing_status"] == "s=-nan"]
        if nan_values:
            ax.scatter([pos] * len(nan_values), nan_values, marker="x", s=58, color="#C9342F", linewidth=1.8, zorder=6)

    ax.axhline(baseline_best, color="#C9453D", linewidth=1.5)
    ax.axvline(1.5, color="#777D81", linestyle=(0, (4, 4)), linewidth=0.9)
    ax.axvline(15.5, color="#777D81", linestyle=(0, (4, 4)), linewidth=0.9)
    ax.set_ylim(ymin, ymax)
    ax.set_xticks(positions)
    ax.set_xticklabels([DISPLAY_NAMES.get(n, n) for n in names], rotation=48, ha="right")
    ax.set_ylabel("Pearson's R", fontsize=11)
    ax.set_title("Figure 3A | Effect of adding family and non-family TFs", fontsize=14, fontweight="bold", pad=28)
    ax.grid(axis="y", color="#D9DDE0", linewidth=0.6, alpha=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    transform = ax.get_xaxis_transform()
    ax.text(1, 1.015, "Baseline", transform=transform, ha="center", fontweight="bold")
    ax.text(8.5, 1.015, "Non-family TFs", transform=transform, ha="center", fontweight="bold")
    ax.text(18.5, 1.015, "ATF/CREB family TFs", transform=transform, ha="center", fontweight="bold")

    for pos, name in zip(positions, names):
        valid_n = len(grouped[name])
        nan_n = sum(row["TF"] == name and row["annealing_status"] == "s=-nan" for row in results)
        label = f"n={valid_n}" + (f" + {nan_n} NaN" if nan_n else "")
        ax.text(pos, ymin + 0.006, label, ha="center", va="bottom", fontsize=6.8, color="#555D62")

    handles = [
        Patch(facecolor=COLORS["baseline"], label="CREB1 self"),
        Patch(facecolor=COLORS["non-family"], label="Non-family TF"),
        Patch(facecolor=COLORS["family"], label="ATF/CREB family TF"),
        plt.Line2D([0], [0], color="#C9453D", linewidth=1.5, label=f"Best baseline R = {baseline_best:.3f}"),
        plt.Line2D([0], [0], marker="x", color="none", markeredgecolor="#C9342F", markeredgewidth=1.8, markersize=7, label="s = -nan run"),
    ]
    ax.legend(handles=handles, frameon=False, loc="lower right", ncol=2)
    fig.subplots_adjust(left=0.065, right=0.99, top=0.88, bottom=0.25)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Reproduce Figure 3A.")
    parser.add_argument("--mmc4", type=Path, default=DEFAULT_MMC4)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--ymin", type=float, default=0.4)
    parser.add_argument("--ymax", type=float, default=0.9)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    results = calculate(args.mmc4)
    metrics = args.output.with_name(f"{args.output.stem}_metrics.csv")
    save_metrics(results, metrics)
    draw(results, args.output, args.ymin, args.ymax, args.dpi, args.show)
    print(f"Saved figure: {args.output}")
    print(f"Saved metrics: {metrics}")


if __name__ == "__main__":
    main()
