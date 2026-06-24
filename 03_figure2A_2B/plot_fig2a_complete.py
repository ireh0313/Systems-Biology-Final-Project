from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd


BASE = Path.home() / "cre_repro"
SCAN_FA = BASE / "work" / "prepared" / "wt_scan_region.fa"
PWM_CALLS = BASE / "work" / "prepared" / "cre_family_pwm_top4_nonoverlap.csv"
CRYPTIC_SUMMARY = BASE / "work" / "prepared" / "cryptic_position_summary.csv"
OUTDIR = BASE / "results" / "fig2"
OUTDIR.mkdir(parents=True, exist_ok=True)

sequence = "".join(
    line.strip() for line in SCAN_FA.read_text().splitlines() if not line.startswith(">")
)
complement = sequence.translate(str.maketrans("ACGT", "TGCA"))

# Physical CRE sites supported by the PWM calls, shown with the concise sequence
# boundaries used for interpretation of the synthetic enhancer.
features = [
    {"label": "CRE1", "start": 11, "end": 19, "color": "#d62728"},
    {"label": "CRE2", "start": 37, "end": 42, "color": "#d62728"},
    {"label": "CRE3", "start": 47, "end": 52, "color": "#d62728"},
    {"label": "cryptic", "start": 63, "end": 67, "color": "#2457d6"},
    {"label": "CRE4", "start": 69, "end": 77, "color": "#d62728"},
]

pwm_calls = pd.read_csv(PWM_CALLS)
cryptic_summary = pd.read_csv(CRYPTIC_SUMMARY)

feature_rows = []
for feature in features:
    row = dict(feature)
    row["sequence"] = sequence[feature["start"] : feature["end"]]
    if feature["label"] == "cryptic":
        evidence = cryptic_summary[
            cryptic_summary["position"].between(feature["start"], feature["end"] - 1)
        ]
        row["evidence"] = (
            "single-hit activity gain + CREB-family PWM score gain; "
            f"max activity gain={evidence['joint_max_activity_gain'].max():.1f}; "
            f"max PWM gain={evidence['joint_max_pwm_gain'].max():.3f}"
        )
    else:
        midpoint = (feature["start"] + feature["end"]) / 2
        calls = pwm_calls.assign(
            midpoint=(pwm_calls["scan_start"] + pwm_calls["scan_end_exclusive"]) / 2
        )
        nearest = calls.iloc[(calls["midpoint"] - midpoint).abs().argsort()[:1]]
        call = nearest.iloc[0]
        row["evidence"] = (
            f"PWM call: {call['TF']} {call['strand']}, score={call['score']:.3f}"
        )
    feature_rows.append(row)

pd.DataFrame(feature_rows).to_csv(
    OUTDIR / "fig2a_complete_feature_coordinates.csv", index=False
)

plt.rcParams.update({"font.family": "DejaVu Sans"})
fig, ax = plt.subplots(figsize=(18, 2.7))

x0 = 3.0
char_step = 1.0
y_top = 0.78
y_bottom = 0.34

ax.text(0.2, y_top, "5′–", ha="left", va="center", fontsize=18)
ax.text(0.2, y_bottom, "3′–", ha="left", va="center", fontsize=18)
ax.text(x0 + len(sequence) * char_step + 0.2, y_top, "–3′", ha="left", va="center", fontsize=18)
ax.text(
    x0 + len(sequence) * char_step + 0.2,
    y_bottom,
    "–5′",
    ha="left",
    va="center",
    fontsize=18,
)

for position, (top_base, bottom_base) in enumerate(zip(sequence, complement)):
    x = x0 + position * char_step
    ax.text(x, y_top, top_base, ha="center", va="center", fontsize=15, family="monospace")
    ax.text(
        x,
        y_bottom,
        bottom_base,
        ha="center",
        va="center",
        fontsize=15,
        family="monospace",
    )

for feature in features:
    start = feature["start"]
    end = feature["end"]
    width = (end - start) * char_step
    left = x0 + start * char_step - 0.5 * char_step
    ax.add_patch(
        Rectangle(
            (left, y_bottom - 0.24),
            width,
            y_top - y_bottom + 0.48,
            fill=False,
            edgecolor=feature["color"],
            linewidth=2.5,
        )
    )
    ax.text(
        left + width / 2,
        y_top + 0.43,
        feature["label"],
        ha="center",
        va="bottom",
        fontsize=18,
        color="black",
        weight="semibold",
    )

ax.set_xlim(-0.5, x0 + len(sequence) + 4)
ax.set_ylim(-0.15, 1.65)
ax.axis("off")
fig.tight_layout()

fig.savefig(OUTDIR / "fig2a_sequence_schematic_complete.png", dpi=300, bbox_inches="tight")
fig.savefig(OUTDIR / "fig2a_sequence_schematic_complete.pdf", bbox_inches="tight")
plt.close(fig)

print("Cryptic region: scan positions 63-66 (63:67), sequence CCCC")
print("Wrote:", OUTDIR / "fig2a_sequence_schematic_complete.png")
print("Wrote:", OUTDIR / "fig2a_sequence_schematic_complete.pdf")
print("Wrote:", OUTDIR / "fig2a_complete_feature_coordinates.csv")
