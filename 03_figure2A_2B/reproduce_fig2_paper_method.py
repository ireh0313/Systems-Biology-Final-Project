from pathlib import Path
import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd


BASE = Path.home() / "cre_repro"
PREP = BASE / "work" / "prepared"
RESULT = BASE / "results" / "fig2" / "paper_method"
RESULT.mkdir(parents=True, exist_ok=True)

sequence = "".join(
    line.strip()
    for line in (PREP / "wt_scan_region.fa").read_text().splitlines()
    if not line.startswith(">")
)
complement = sequence.translate(str.maketrans("ACGT", "TGCA"))

# CRE sites are known properties of the designed synthetic enhancer.
canonical = pd.DataFrame(
    [
        ("CRE1", 11, 19, "TGACGTCA", "full-length CRE"),
        ("CRE2", 37, 42, "CGTCA", "half-length CRE"),
        ("CRE3", 47, 52, "TGACG", "reverse-strand half-length CRE (CGTCA)"),
        ("CRE4", 69, 77, "TGACGTCA", "full-length CRE"),
    ],
    columns=["label", "start", "end_exclusive", "display_sequence", "definition"],
)
canonical["observed_sequence"] = canonical.apply(
    lambda row: sequence[row.start : row.end_exclusive], axis=1
)
canonical.to_csv(RESULT / "canonical_cre_coordinates.csv", index=False)

# Reproduce the paper's threshold-based use of PWM sites. The raw PWM scan was
# already performed at 1-bp shifts on both strands; here every hit is filtered
# using the corresponding fitted 7TF threshold instead of selecting top four.
hits = pd.read_csv(PREP / "cre_family_pwm_all_hits.csv")
xml_path = BASE.parent / "transcpp-master" / "CRE_7TF_self8_singlehit_fullanneal.xml"
root = ET.parse(xml_path).getroot()
thresholds = {}
for tf_node in root.findall(".//TF"):
    threshold_node = tf_node.find("threshold")
    if threshold_node is not None:
        thresholds[tf_node.attrib["name"]] = float(threshold_node.attrib["value"])

hits["threshold"] = hits["TF"].map(thresholds)
hits["passes_threshold"] = hits["score"] >= hits["threshold"]
passing = hits[hits["passes_threshold"]].copy()
passing.to_csv(RESULT / "threshold_passing_pwm_sites.csv", index=False)

validation_rows = []
for _, site in canonical.iterrows():
    overlaps = passing[
        (passing["scan_start"] < site.end_exclusive)
        & (passing["scan_end_exclusive"] > site.start)
    ].sort_values("score", ascending=False)
    validation_rows.append(
        {
            "canonical_site": site.label,
            "canonical_start": site.start,
            "canonical_end_exclusive": site.end_exclusive,
            "canonical_sequence": site.observed_sequence,
            "n_threshold_passing_overlaps": len(overlaps),
            "best_overlap_TF": overlaps.iloc[0]["TF"] if len(overlaps) else "",
            "best_overlap_strand": overlaps.iloc[0]["strand"] if len(overlaps) else "",
            "best_overlap_start": overlaps.iloc[0]["scan_start"] if len(overlaps) else "",
            "best_overlap_end_exclusive": (
                overlaps.iloc[0]["scan_end_exclusive"] if len(overlaps) else ""
            ),
            "best_overlap_score": overlaps.iloc[0]["score"] if len(overlaps) else "",
            "best_overlap_threshold": (
                overlaps.iloc[0]["threshold"] if len(overlaps) else ""
            ),
        }
    )
pd.DataFrame(validation_rows).to_csv(
    RESULT / "canonical_cre_pwm_threshold_validation.csv", index=False
)

# Paper-like cryptic discovery: use unexpected measured MPRA activity outside
# known CRE sites, without requiring a PWM gain.
data = pd.read_csv(PREP / "single_hit_delta.csv")
mutants = data[data["kind"] == "single-hit"].copy()
mutants["position"] = mutants["position"].astype(int)
mutants["delta_activity_log2"] = mutants["relative_expression"].apply(
    lambda value: __import__("math").log2(value)
)
position_summary = (
    mutants.groupby("position")
    .agg(
        mean_delta_log2=("delta_activity_log2", "mean"),
        max_delta_log2=("delta_activity_log2", "max"),
        min_delta_log2=("delta_activity_log2", "min"),
        sd_delta_log2=("delta_activity_log2", "std"),
        n_positive_log2=(
            "delta_activity_log2",
            lambda values: int((values > 0).sum()),
        ),
    )
    .reset_index()
)
position_summary["in_canonical_cre"] = position_summary["position"].apply(
    lambda position: any(
        row.start <= position < row.end_exclusive for _, row in canonical.iterrows()
    )
)
position_summary["unexpected_positive"] = (
    (~position_summary["in_canonical_cre"])
    & (position_summary["n_positive_log2"] == 3)
    & (position_summary["mean_delta_log2"] >= 0.1)
    & (position_summary["sd_delta_log2"] <= 0.15)
)
position_summary.to_csv(RESULT / "mpra_position_summary.csv", index=False)

candidate_positions = position_summary[
    position_summary["unexpected_positive"]
]["position"].tolist()
runs = []
for position in candidate_positions:
    if not runs or position != runs[-1][-1] + 1:
        runs.append([position])
    else:
        runs[-1].append(position)
cryptic_run = max(runs, key=lambda run: (len(run), sum(run)))
cryptic_start = cryptic_run[0]
cryptic_end = cryptic_run[-1] + 1
cryptic = pd.DataFrame(
    [
        {
            "label": "cryptic candidate",
            "start": cryptic_start,
            "end_exclusive": cryptic_end,
            "sequence": sequence[cryptic_start:cryptic_end],
            "rule": (
                "longest contiguous non-CRE run where all 3 substitutions "
                "have positive log2(mutant/WT), mean log2 delta >= 0.1, "
                "and substitution-effect SD <= 0.15"
            ),
        }
    ]
)
cryptic.to_csv(RESULT / "cryptic_candidate_mpra_only.csv", index=False)


def add_spans(ax):
    for _, site in canonical.iterrows():
        ax.axvspan(site.start - 0.5, site.end_exclusive - 0.5, color="#ef9a9a", alpha=0.45)
    ax.axvspan(cryptic_start - 0.5, cryptic_end - 0.5, color="#9fa8ff", alpha=0.55)


# Paper-like Figure 2B: four substitution-specific MPRA bar panels.
# The manuscript displays delta activity as log2(mutant activity / WT activity).
bases = ["A", "C", "G", "T"]
pivot = mutants.pivot(
    index="position", columns="mut_base", values="delta_activity_log2"
)
fig, axes = plt.subplots(4, 1, figsize=(13, 9), sharex=True, sharey=True)
for ax, base in zip(axes, bases):
    add_spans(ax)
    values = pivot[base].dropna()
    ax.bar(values.index, values.values, width=0.75, color="#123cff")
    ax.axhline(0, color="black", linewidth=1.5)
    ax.text(0.02, 0.82, f"→{base}", transform=ax.transAxes, fontsize=18)
    ax.set_ylim(-2, 2)
    ax.set_yticks([-2, -1, 0, 1, 2])
    ax.grid(axis="both", color="#e8e8e8", linewidth=0.8)
    ax.set_axisbelow(True)
axes[0].set_title("Measured mRNA", fontsize=20, color="black", pad=12)
axes[-1].set_xlabel("position", fontsize=16)
axes[-1].tick_params(axis="x", labelbottom=False)
fig.supylabel("Δactivity (log2)", fontsize=16)
fig.tight_layout()
fig.savefig(RESULT / "fig2b_single_hit_mpra_paper_method.png", dpi=300)
fig.savefig(RESULT / "fig2b_single_hit_mpra_paper_method.pdf")
plt.close(fig)

# Paper-like Figure 2A using known canonical sites and independently detected
# MPRA-only cryptic candidate.
features = [
    *[
        {
            "label": row.label,
            "start": int(row.start),
            "end": int(row.end_exclusive),
            "color": "#d62728",
        }
        for _, row in canonical.iterrows()
    ],
    {
        "label": "cryptic candidate",
        "start": cryptic_start,
        "end": cryptic_end,
        "color": "#2457d6",
    },
]
features.sort(key=lambda feature: feature["start"])

fig, ax = plt.subplots(figsize=(18, 2.7))
x0, y_top, y_bottom = 3.0, 0.78, 0.34
ax.text(0.2, y_top, "5′–", ha="left", va="center", fontsize=18)
ax.text(0.2, y_bottom, "3′–", ha="left", va="center", fontsize=18)
ax.text(x0 + len(sequence) + 0.2, y_top, "–3′", ha="left", va="center", fontsize=18)
ax.text(x0 + len(sequence) + 0.2, y_bottom, "–5′", ha="left", va="center", fontsize=18)
for position, (top_base, bottom_base) in enumerate(zip(sequence, complement)):
    x = x0 + position
    ax.text(x, y_top, top_base, ha="center", va="center", fontsize=15, family="monospace")
    ax.text(x, y_bottom, bottom_base, ha="center", va="center", fontsize=15, family="monospace")
for feature in features:
    left = x0 + feature["start"] - 0.5
    width = feature["end"] - feature["start"]
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
        fontsize=16,
        color="black",
        weight="semibold",
    )
ax.set_xlim(-0.5, x0 + len(sequence) + 4)
ax.set_ylim(-0.15, 1.65)
ax.axis("off")
fig.tight_layout()
fig.savefig(RESULT / "fig2a_sequence_paper_method.png", dpi=300, bbox_inches="tight")
fig.savefig(RESULT / "fig2a_sequence_paper_method.pdf", bbox_inches="tight")
plt.close(fig)

print(f"Threshold-passing PWM sites: {len(passing)}")
print(f"Cryptic candidate: {cryptic_start}:{cryptic_end} {sequence[cryptic_start:cryptic_end]}")
print(f"Output directory: {RESULT}")
