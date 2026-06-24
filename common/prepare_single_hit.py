from pathlib import Path
import re
import pandas as pd

base = Path.home() / "cre_repro"
inp = base / "work" / "csv" / "seqs__single-hit.csv"
outdir = base / "work" / "prepared"
outdir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(inp)
df["expression level"] = pd.to_numeric(df["expression level"])
df["sequence_length"] = df["sequence"].str.len()

wt = df.iloc[0]
wt_name = wt["name"]
wt_seq = wt["sequence"]
wt_expr = float(wt["expression level"])

records = []
records.append({
    "name": wt_name,
    "kind": "WT",
    "position": None,
    "mut_base": None,
    "wt_base": None,
    "sequence": wt_seq,
    "expression": wt_expr,
    "delta_expression": 0.0,
    "relative_expression": 1.0,
    "sequence_length": len(wt_seq),
})

pat = re.compile(r"scanmut_single_pos_(\d+)_([ACGT])$")
for _, row in df.iloc[1:].iterrows():
    name = row["name"]
    m = pat.match(name)
    if not m:
        raise ValueError(f"Unexpected single-hit name: {name}")

    pos = int(m.group(1))
    mut_base = m.group(2)
    seq = row["sequence"]
    expr = float(row["expression level"])

    records.append({
        "name": name,
        "kind": "single-hit",
        "position": pos,
        "mut_base": mut_base,
        "wt_base": wt_seq[20 + pos] if 20 + pos < len(wt_seq) else "",
        "sequence": seq,
        "expression": expr,
        "delta_expression": expr - wt_expr,
        "relative_expression": expr / wt_expr,
        "sequence_length": len(seq),
    })

out = pd.DataFrame(records)
out.to_csv(outdir / "single_hit_delta.csv", index=False)

summary = {
    "wt_name": wt_name,
    "wt_expression": wt_expr,
    "n_total": len(out),
    "n_single_hit": int((out["kind"] == "single-hit").sum()),
    "sequence_lengths": sorted(out["sequence_length"].unique().tolist()),
    "positions_min": int(out["position"].dropna().min()),
    "positions_max": int(out["position"].dropna().max()),
}
pd.Series(summary).to_csv(outdir / "single_hit_summary.csv", header=False)

print("Wrote:", outdir / "single_hit_delta.csv")
print("Wrote:", outdir / "single_hit_summary.csv")
print(summary)
print(out.head())
