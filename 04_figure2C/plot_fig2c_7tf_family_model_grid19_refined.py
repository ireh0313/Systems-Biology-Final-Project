from pathlib import Path
import re
import subprocess

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd


HOME = Path.home()
BASE = HOME / "cre_repro"
XML = BASE / "results/aws_7TF_full/CRE_7TF_full_seed870001.xml"
UNFOLD = HOME / "transcpp-master/unfold"
MEASURED = BASE / "work/prepared/single_hit_delta.csv"
OUT = BASE / "results/fig2/paper_method"

canonical = [(11, 19), (37, 42), (47, 52), (69, 77)]
cryptic = (63, 67)
name_pattern = re.compile(r"scanmut_single_pos_(\d+)_([ACGT])$")


def read_predicted_rates():
    text = subprocess.check_output(
        [str(UNFOLD), "-i", str(XML), "--rate", "--invert"],
        text=True,
    )
    rows = []
    for line in text.splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[0] != "id":
            rows.append((parts[0], float(parts[1])))
    return pd.DataFrame(rows, columns=["name", "predicted_rate"])


predicted = read_predicted_rates()
predicted_wt = predicted.loc[
    predicted["name"] == "synCRE_Promega_0", "predicted_rate"
].iloc[0]
predicted["predicted_delta_log2"] = np.log2(predicted["predicted_rate"] / predicted_wt)

measured = pd.read_csv(MEASURED)
measured = measured[measured["kind"] == "single-hit"].copy()
measured["measured_delta_log2"] = np.log2(measured["relative_expression"])

combined = measured.merge(
    predicted[["name", "predicted_rate", "predicted_delta_log2"]],
    on="name",
    how="inner",
)
combined["position"] = combined["position"].astype(int)


def add_spans(ax):
    for start, end in canonical:
        ax.axvspan(start - 0.5, end - 0.5, color="#ef9a9a", alpha=0.45)
    ax.axvspan(cryptic[0] - 0.5, cryptic[1] - 0.5, color="#9fa8ff", alpha=0.55)


# Keep 19 total vertical grid cells while leaving one empty cell at each side.
data_left, data_right = -0.5, 86.5
cell_width = (data_right - data_left) / 17
plot_left = data_left - cell_width
plot_right = data_right + cell_width
grid_boundaries = np.linspace(plot_left, plot_right, 20)

fig, axes = plt.subplots(4, 1, figsize=(8.4, 10.0), sharex=True, sharey=True)
metrics = []

for ax, base in zip(axes, ["A", "C", "G", "T"]):
    subset = combined[combined["mut_base"] == base].sort_values("position")
    observed = subset["measured_delta_log2"].to_numpy()
    modeled = subset["predicted_delta_log2"].to_numpy()

    r = float(np.corrcoef(observed, modeled)[0, 1])
    r_squared = r**2
    rmse = float(np.sqrt(np.mean((observed - modeled) ** 2)))
    metrics.append((base, len(subset), r, r_squared, rmse))

    add_spans(ax)
    ax.bar(subset["position"], modeled, width=0.96, color="#f00000")
    ax.hlines(0, data_left, data_right, color="black", linewidth=1.8)
    ax.text(plot_left + 4.0, 1.52, f"\u2192{base}", fontsize=18, va="center")
    ax.text(
        data_right - 1.0,
        1.52,
        rf"R= {r:.2f}   R$^2$= {r_squared:.2f}   rmse= {rmse:.2f}",
        ha="right",
        va="center",
        fontsize=17,
        color="black",
    )

    ax.set_xlim(plot_left, plot_right)
    ax.set_ylim(-2.25, 2.25)
    ax.set_yticks([-2, -1, 0, 1, 2])
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax.set_xticks(np.arange(10, 81, 10))
    ax.set_xticks(grid_boundaries, minor=True)
    ax.tick_params(axis="x", labelbottom=False)
    ax.tick_params(axis="x", which="minor", length=0)

    ax.grid(which="major", axis="y", color="#e5e5e5", linewidth=0.8)
    ax.grid(which="minor", axis="y", color="#e5e5e5", linewidth=0.8)
    ax.grid(which="minor", axis="x", color="#e5e5e5", linewidth=0.8)
    ax.set_axisbelow(True)

    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
        spine.set_color("black")

axes[-1].tick_params(axis="x", labelbottom=True, labelsize=12)
axes[0].set_title("7 CREB family model", fontsize=22, color="black", pad=12)
axes[-1].set_xlabel("position", fontsize=17)
fig.supylabel(r"predicted $\Delta$activity (log2)", fontsize=17, x=0.01)
fig.subplots_adjust(left=0.12, right=0.985, top=0.94, bottom=0.07, hspace=0.13)

png = OUT / "fig2c_7tf_family_model_grid19_refined.png"
pdf = OUT / "fig2c_7tf_family_model_grid19_refined.pdf"
csv = OUT / "fig2c_7tf_family_model_measured_vs_predicted.csv"
metric_file = OUT / "fig2c_7tf_family_model_metrics.csv"

fig.savefig(png, dpi=300)
fig.savefig(pdf)
plt.close(fig)

combined.to_csv(csv, index=False)
pd.DataFrame(
    metrics, columns=["mut_base", "n", "pearson_r", "r_squared", "rmse_log2"]
).to_csv(metric_file, index=False)

print(png)
print(pdf)
print(metric_file)
