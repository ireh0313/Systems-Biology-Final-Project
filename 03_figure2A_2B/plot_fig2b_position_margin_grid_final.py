from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd


BASE = Path.home() / "cre_repro"
PREP = BASE / "work" / "prepared"
OUT = BASE / "results" / "fig2" / "paper_method"

canonical = [(11, 19), (37, 42), (47, 52), (69, 77)]
cryptic = (63, 67)

data = pd.read_csv(PREP / "single_hit_delta.csv")
mutants = data[data["kind"] == "single-hit"].copy()
mutants["position"] = mutants["position"].astype(int)
mutants["delta_activity_log2"] = np.log2(mutants["relative_expression"])
pivot = mutants.pivot(index="position", columns="mut_base", values="delta_activity_log2")


def add_spans(ax):
    for start, end in canonical:
        ax.axvspan(start - 0.5, end - 0.5, color="#ef9a9a", alpha=0.45)
    ax.axvspan(cryptic[0] - 0.5, cryptic[1] - 0.5, color="#9fa8ff", alpha=0.55)


def add_paper_annotations(ax, base):
    arrow = dict(arrowstyle="-|>", mutation_scale=17, linewidth=3.0)
    ax.annotate("", xy=(14, -0.78), xytext=(10, -1.22),
                arrowprops={**arrow, "color": "#00a000"})

    if base != "T":
        ax.annotate("", xy=(67, -0.70), xytext=(63, -1.12),
                    arrowprops={**arrow, "color": "#f00000"})
        ax.annotate("", xy=(78, -0.72), xytext=(82, -1.15),
                    arrowprops={**arrow, "color": "#f00000"})
    else:
        ax.annotate("", xy=(76.5, -0.96), xytext=(80.5, -1.39),
                    arrowprops={**arrow, "color": "#f00000"})

    if base != "C":
        ax.text(64.5, 0.82, "*", ha="center", va="center", fontsize=30, color="black")


# Keep 19 total grid cells. The data region occupies the inner 17 cells,
# leaving exactly one empty grid cell at both the left and right edges.
data_left, data_right = -0.5, 86.5
cell_width = (data_right - data_left) / 17
plot_left = data_left - cell_width
plot_right = data_right + cell_width
grid_boundaries = np.linspace(plot_left, plot_right, 20)

fig, axes = plt.subplots(4, 1, figsize=(7.2, 10.0), sharex=True, sharey=True)

for ax, base in zip(axes, ["A", "C", "G", "T"]):
    add_spans(ax)
    values = pivot[base].dropna()
    ax.bar(values.index, values.values, width=0.96, color="#123cff")
    # The black zero baseline represents only the measured sequence region;
    # leave the two outer grid cells as true visual margins.
    ax.hlines(0, data_left, data_right, color="black", linewidth=1.7)
    ax.text(0.035, 0.80, f"\u2192{base}", transform=ax.transAxes, fontsize=18)

    ax.set_xlim(plot_left, plot_right)
    ax.set_ylim(-2.25, 2.25)
    ax.set_yticks([-2, -1, 0, 1, 2])
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))

    # Major x ticks carry the requested numeric labels; minor ticks draw the
    # 20 boundaries that form exactly 19 grid cells.
    ax.set_xticks(np.arange(0, 81, 10))
    ax.set_xticks(grid_boundaries, minor=True)
    ax.tick_params(axis="x", which="major", labelbottom=False, length=4)
    ax.tick_params(axis="x", which="minor", length=0)

    ax.grid(which="minor", axis="y", color="#e5e5e5", linewidth=0.8)
    ax.grid(which="major", axis="y", color="#e5e5e5", linewidth=0.8)
    ax.grid(which="minor", axis="x", color="#e5e5e5", linewidth=0.8)
    ax.set_axisbelow(True)
    add_paper_annotations(ax, base)

axes[-1].tick_params(axis="x", which="major", labelbottom=True, labelsize=11)
axes[0].set_title("Measured mRNA", fontsize=21, color="black", pad=12)
axes[-1].set_xlabel("position", fontsize=17)
fig.supylabel(r"$\Delta$activity (log2)", fontsize=17, x=0.01)
fig.subplots_adjust(left=0.13, right=0.985, top=0.94, bottom=0.07, hspace=0.13)

png = OUT / "fig2b_single_hit_mpra_position_margin_grid_final.png"
pdf = OUT / "fig2b_single_hit_mpra_position_margin_grid_final.pdf"
fig.savefig(png, dpi=300)
fig.savefig(pdf)
plt.close(fig)
print(png)
print(pdf)
