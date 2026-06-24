from pathlib import Path
import numpy as np
import pandas as pd


BASE = Path.home() / "cre_repro"
DATA = BASE / "work" / "prepared" / "single_hit_delta.csv"
OUT = BASE / "results" / "fig2" / "paper_method"
OUT.mkdir(parents=True, exist_ok=True)

canonical = [(11, 19), (37, 42), (47, 52), (69, 77)]
data = pd.read_csv(DATA)
mut = data[data["kind"] == "single-hit"].copy()
mut["position"] = mut["position"].astype(int)
mut["delta_log2"] = np.log2(mut["relative_expression"])

summary = (
    mut.groupby(["position", "wt_base"])
    .agg(
        mean_log2=("delta_log2", "mean"),
        min_log2=("delta_log2", "min"),
        max_log2=("delta_log2", "max"),
        sd_log2=("delta_log2", "std"),
        n_positive=("delta_log2", lambda values: int((values > 0).sum())),
    )
    .reset_index()
)
summary["in_canonical_cre"] = summary["position"].apply(
    lambda position: any(start <= position < end for start, end in canonical)
)


def contiguous_runs(positions):
    runs = []
    for position in sorted(positions):
        if not runs or position != runs[-1][-1] + 1:
            runs.append([position])
        else:
            runs[-1].append(position)
    return runs


def sequence_for(run):
    sequence = "".join(
        line.strip()
        for line in (BASE / "work" / "prepared" / "wt_scan_region.fa").read_text().splitlines()
        if not line.startswith(">")
    )
    return sequence[run[0] : run[-1] + 1]


rows = []
mean_thresholds = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4]
sd_thresholds = [None, 0.25, 0.2, 0.15, 0.1, 0.05]

for mean_threshold in mean_thresholds:
    for sd_threshold in sd_thresholds:
        selected = summary[
            (~summary["in_canonical_cre"])
            & (summary["n_positive"] == 3)
            & (summary["mean_log2"] >= mean_threshold)
        ]
        if sd_threshold is not None:
            selected = selected[selected["sd_log2"] <= sd_threshold]

        runs = contiguous_runs(selected["position"].tolist())
        longest = max(runs, key=lambda run: (len(run), sum(run))) if runs else []
        rows.append(
            {
                "mean_log2_min": mean_threshold,
                "sd_log2_max": "none" if sd_threshold is None else sd_threshold,
                "selected_positions": ";".join(map(str, selected["position"].tolist())),
                "all_runs": ";".join(f"{run[0]}:{run[-1] + 1}" for run in runs),
                "longest_start": longest[0] if longest else "",
                "longest_end_exclusive": longest[-1] + 1 if longest else "",
                "longest_sequence": sequence_for(longest) if longest else "",
            }
        )

result = pd.DataFrame(rows)
result.to_csv(OUT / "cryptic_threshold_sensitivity.csv", index=False)

print("All mean-only criteria:")
print(result[result["sd_log2_max"] == "none"].to_string(index=False))
print("\nCriteria yielding longest sequence CCCC:")
print(result[result["longest_sequence"] == "CCCC"].to_string(index=False))
