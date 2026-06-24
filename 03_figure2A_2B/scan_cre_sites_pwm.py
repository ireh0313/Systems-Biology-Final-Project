from pathlib import Path
import pandas as pd

base = Path.home() / "cre_repro"
pwms_path = base / "work" / "csv" / "pwms.csv"
single_path = base / "work" / "prepared" / "single_hit_delta.csv"
outdir = base / "work" / "prepared"
outdir.mkdir(parents=True, exist_ok=True)

tf_keep = ["CREB1", "CREB3", "CREB5", "CREM", "ATF1", "ATF4", "ATF7"]
base_index = {"A": 0, "C": 1, "G": 2, "T": 3}
comp = str.maketrans("ACGT", "TGCA")

def revcomp(s):
    return s.translate(comp)[::-1]

# WT scan region
df = pd.read_csv(single_path)
wt_seq = df[df["kind"] == "WT"].iloc[0]["sequence"]
scan_start = 20
scan_end = 107
scan_seq = wt_seq[scan_start:scan_end]

# Parse PWM/PSSM blocks
raw = pd.read_csv(pwms_path, header=None, dtype=str).fillna("")
pwms = {}

i = 0
while i < len(raw):
    tf = raw.iat[i, 0]
    letter = raw.iat[i, 1]
    if tf in tf_keep and letter == "A":
        rows = raw.iloc[i:i+4]
        letters = list(rows.iloc[:, 1])
        if letters != ["A", "C", "G", "T"]:
            raise ValueError(f"Unexpected PWM rows for {tf}: {letters}")

        mat = rows.iloc[:, 2:].replace("", pd.NA).dropna(axis=1, how="all").astype(float)
        pwms[tf] = mat.to_numpy()
        i += 4
    else:
        i += 1

# Score all windows on both strands
hits = []
for tf, mat in pwms.items():
    L = mat.shape[1]
    for start in range(0, len(scan_seq) - L + 1):
        window = scan_seq[start:start+L]

        score_f = sum(mat[base_index[b], j] for j, b in enumerate(window))
        hits.append({
            "TF": tf,
            "strand": "+",
            "scan_start": start,
            "scan_end_exclusive": start + L,
            "construct_start": scan_start + start,
            "construct_end_exclusive": scan_start + start + L,
            "width": L,
            "sequence": window,
            "score": score_f,
        })

        rc = revcomp(window)
        score_r = sum(mat[base_index[b], j] for j, b in enumerate(rc))
        hits.append({
            "TF": tf,
            "strand": "-",
            "scan_start": start,
            "scan_end_exclusive": start + L,
            "construct_start": scan_start + start,
            "construct_end_exclusive": scan_start + start + L,
            "width": L,
            "sequence": window,
            "score": score_r,
        })

hits = pd.DataFrame(hits)
hits = hits.sort_values("score", ascending=False)
hits.to_csv(outdir / "cre_family_pwm_all_hits.csv", index=False)

# Pick non-overlapping top CRE-like clusters.
# This prevents the same physical site from being counted many times by different TFs.
chosen = []
for _, row in hits.iterrows():
    start = int(row["scan_start"])
    end = int(row["scan_end_exclusive"])

    overlaps = False
    for c in chosen:
        if not (end <= c["scan_start"] or start >= c["scan_end_exclusive"]):
            overlaps = True
            break

    if not overlaps:
        chosen.append(row.to_dict())

    if len(chosen) == 4:
        break

chosen_df = pd.DataFrame(chosen)
chosen_df.insert(0, "CRE_call", [f"CRE-like_{i}" for i in range(1, len(chosen_df)+1)])
chosen_df.to_csv(outdir / "cre_family_pwm_top4_nonoverlap.csv", index=False)

print("WT scan sequence:")
print(scan_seq)
print()
print("Top 20 PWM hits:")
print(hits.head(20)[["TF", "strand", "scan_start", "scan_end_exclusive", "sequence", "score"]].to_string(index=False))
print()
print("Top 4 non-overlapping CRE-like calls:")
print(chosen_df[["CRE_call", "TF", "strand", "scan_start", "scan_end_exclusive", "sequence", "score"]].to_string(index=False))
print()
print("Wrote:")
print(outdir / "cre_family_pwm_all_hits.csv")
print(outdir / "cre_family_pwm_top4_nonoverlap.csv")
