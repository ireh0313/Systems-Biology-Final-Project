from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HOME = Path.home()
RATE_FILE = HOME / "CRE_1TF_selfoff_singlehit_fullanneal_CREB1_scramble_rate_invert.txt"
OUT = HOME / "cre_repro/results/fig2/paper_method"

rows = []
for line in RATE_FILE.read_text().splitlines():
    parts = line.split()
    if len(parts) != 2 or parts[0] == "id":
        continue
    rows.append((parts[0], float(parts[1])))

rates = pd.DataFrame(rows, columns=["name", "predicted_rate"])
wt_rate = rates.loc[rates["name"] == "synCRE_Promega_0", "predicted_rate"].iloc[0]

pattern = re.compile(r"scanmut_single_pos_(\d+)_([ACGT])$")
records = []
for row in rates.itertuples(index=False):
    match = pattern.match(row.name)
    if match:
        records.append(
            {
                "position": int(match.group(1)),
                "mut_base": match.group(2),
                "delta_activity_log2": np.log2(row.predicted_rate / wt_rate),
            }
        )

data = pd.DataFrame(records)
pivot = data.pivot(index="position", columns="mut_base", values="delta_activity_log2")

fig, axes = plt.subplots(4, 1, figsize=(7.2, 7.0), sharex=True, sharey=True)

for ax, base in zip(axes, ["A", "C", "G", "T"]):
    values = pivot[base].dropna()
    ax.bar(values.index, values.values, width=0.96, color="#f00000")
    ax.axhline(0, color="black", linewidth=1.8)

    ax.set_xlim(-0.5, 86.5)
    ax.set_ylim(-2.25, 0.55)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.tick_params(
        axis="both",
        which="both",
        bottom=False,
        top=False,
        left=False,
        right=False,
        labelbottom=False,
        labelleft=False,
    )

    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(1.6)

axes[0].text(0.02, 0.08, "B.M", transform=axes[0].transAxes, fontsize=14, color="#555555")
fig.subplots_adjust(left=0.02, right=0.995, top=0.995, bottom=0.01, hspace=0.14)

png = OUT / "fig2c_creb1_baseline_model_clean_inset.png"
pdf = OUT / "fig2c_creb1_baseline_model_clean_inset.pdf"
fig.savefig(png, dpi=300)
fig.savefig(pdf)
plt.close(fig)

summary = OUT / "fig2c_creb1_baseline_model_summary.txt"
summary.write_text(
    f"source={RATE_FILE}\n"
    f"predicted_wt={wt_rate}\n"
    f"n_variants={len(data)}\n"
    f"min_log2_delta={data.delta_activity_log2.min()}\n"
    f"max_log2_delta={data.delta_activity_log2.max()}\n"
)

print(png)
print(pdf)
print(summary)
