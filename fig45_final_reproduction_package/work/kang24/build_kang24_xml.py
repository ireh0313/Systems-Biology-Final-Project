#!/usr/bin/env python3
"""Build transcpp XML inputs from Kang and Kim 2024 supplemental tables."""

from __future__ import annotations

import argparse
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd


BASES = ("A", "C", "G", "T")
SUPPLEMENTARY = Path(__file__).resolve().parent / "supplementary"


def valid_xml_name(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]", "_", str(name))
    if not re.match(r"^[A-Za-z_]", cleaned):
        cleaned = f"g_{cleaned}"
    return cleaned


def bool_text(value: object) -> str:
    if isinstance(value, str):
        return value.lower()
    return "true" if bool(value) else "false"


def fmt(value: object) -> str:
    if isinstance(value, str):
        return value
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return f"{float(value):.12g}"


def param(parent: ET.Element, tag: str, value: object, lim_low: object, lim_high: object, anneal: object, move: str | None = None) -> ET.Element:
    attrs = {
        "value": fmt(value),
        "lim_low": fmt(lim_low),
        "lim_high": fmt(lim_high),
        "anneal": bool_text(anneal),
    }
    if move:
        attrs["move"] = move
    return ET.SubElement(parent, tag, attrs)


def read_model_parameters(path: Path, sheet: str, model_column: int | None) -> tuple[dict, dict[str, dict]]:
    df = pd.read_excel(path, sheet_name=sheet, header=None)
    score_rows = df.index[df.iloc[:, 1].astype(str).str.contains("score", case=False, na=False)]
    score_row = score_rows[-1]
    candidate_cols = [
        c
        for c in range(2, min(10, df.shape[1]))
        if pd.notna(df.iat[score_row, c])
    ]
    if model_column is None:
        col = min(candidate_cols, key=lambda c: float(df.iat[score_row, c]))
    else:
        col = model_column

    limits_low_col = 10
    limits_high_col = 11
    tweak_col = 12
    promoter = {}
    tfs: dict[str, dict] = {}

    for _, row in df.iloc[10 : score_row + 1].iterrows():
        tf = row.iloc[0]
        field = row.iloc[1]
        if pd.isna(field) or str(field).startswith("score"):
            continue
        field = str(field)
        entry = {
            "value": row.iloc[col],
            "lim_low": row.iloc[limits_low_col],
            "lim_high": row.iloc[limits_high_col],
            "anneal": row.iloc[tweak_col],
        }
        if pd.isna(tf):
            promoter[field] = entry
        else:
            tfs.setdefault(str(tf), {})[field] = entry

    promoter["model_num"] = df.iat[9, col]
    promoter["score_sse"] = df.iat[score_row, col]
    promoter["sheet"] = sheet
    promoter["source_column"] = col
    return promoter, tfs


def read_pwms(path: Path) -> dict[str, list[list[float]]]:
    df = pd.read_excel(path, sheet_name=0, header=None)
    pwms: dict[str, list[list[float]]] = {}
    for row in range(df.shape[0] - 3):
        tf = df.iat[row, 0]
        base = df.iat[row, 1]
        if pd.isna(tf) or base != "A":
            continue
        base_rows = {}
        for offset, expected_base in enumerate(BASES):
            if df.iat[row + offset, 1] != expected_base:
                raise ValueError(f"Unexpected PWM block for {tf}: missing {expected_base}")
            values = [
                float(v)
                for v in df.iloc[row + offset, 2:].tolist()
                if pd.notna(v)
            ]
            base_rows[expected_base] = values
        length = len(base_rows["A"])
        if any(len(base_rows[b]) != length for b in BASES):
            raise ValueError(f"Unequal PWM row lengths for {tf}")
        pwms[str(tf)] = [[base_rows[b][i] for b in BASES] for i in range(length)]
    return pwms


def read_sequences(path: Path, sheet: str, limit: int | None) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name=sheet)
    required = {"name", "sequence", "expression level"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{sheet} is missing columns: {sorted(missing)}")
    df = df.dropna(subset=["name", "sequence", "expression level"]).copy()
    df["xml_name"] = [valid_xml_name(v) for v in df["name"]]
    if df["xml_name"].duplicated().any():
        counts = {}
        names = []
        for name in df["xml_name"]:
            counts[name] = counts.get(name, 0) + 1
            names.append(name if counts[name] == 1 else f"{name}_{counts[name]}")
        df["xml_name"] = names
    if limit:
        df = df.iloc[:limit].copy()
    return df


def add_annealer(root: ET.Element, smoke: bool, init_loop: int | None) -> None:
    if smoke:
        ET.SubElement(root, "annealer_input", {"init_T": "100000", "lambda": "0.0001", "init_loop": str(init_loop or 10)})
        ET.SubElement(root, "move", {"interval": "10", "gain": "3"})
        ET.SubElement(root, "count_criterion", {"freeze_crit": "1", "freeze_cnt": "1"})
        ET.SubElement(root, "mix", {"interval": "10", "adaptcoef": "10"})
        ET.SubElement(root, "lam", {"tau": "10", "memLength_mean": ".200", "memLength_sd": "10", "criterion": "1", "freeze_cnt": "1"})
    else:
        ET.SubElement(root, "annealer_input", {"init_T": "100000", "lambda": "0.0001", "init_loop": str(init_loop or 100000)})
        ET.SubElement(root, "move", {"interval": "100", "gain": "3"})
        ET.SubElement(root, "count_criterion", {"freeze_crit": "1", "freeze_cnt": "5"})
        ET.SubElement(root, "mix", {"interval": "100", "adaptcoef": "10"})
        ET.SubElement(root, "lam", {"tau": "100", "memLength_mean": ".200", "memLength_sd": "100", "criterion": "1", "freeze_cnt": "5"})


def build_xml(
    promoter: dict,
    tfs: dict[str, dict],
    pwms: dict[str, list[list[float]]],
    seqs: pd.DataFrame,
    smoke: bool,
    fixed: bool,
    seed: str,
    init_loop: int | None,
    competition: bool,
    self_competition: bool,
) -> ET.ElementTree:
    root = ET.Element("Root")
    add_annealer(root, smoke, init_loop)

    mode = ET.SubElement(root, "Mode")
    ET.SubElement(mode, "Schedule", {"value": "0"})
    ET.SubElement(mode, "Verbose", {"value": "1" if smoke else "2"})
    ET.SubElement(mode, "Profiling", {"value": "false"})
    ET.SubElement(mode, "ScoreFunction", {"value": "sse"})
    ET.SubElement(mode, "PThresh", {"value": "0"})
    ET.SubElement(mode, "ScaleData", {"value": "false", "type": "area", "scale_to": "200"})
    ET.SubElement(mode, "PerGene", {"value": "false"})
    ET.SubElement(mode, "PerNuc", {"value": "false"})
    ET.SubElement(mode, "MinData", {"value": "0"})
    ET.SubElement(mode, "Competition", {"value": bool_text(competition), "window": "500", "shift": "50"})
    ET.SubElement(mode, "NumThreads", {"value": "4"})
    ET.SubElement(mode, "SelfCompetition", {"value": bool_text(self_competition)})
    ET.SubElement(mode, "Precision", {"value": "16"})
    ET.SubElement(mode, "Seed", {"value": seed})
    ET.SubElement(mode, "GCcontent", {"value": "0.5"})

    inp = ET.SubElement(root, "Input")
    distances = ET.SubElement(inp, "Distances")
    quenching = ET.SubElement(distances, "Distance", {"name": "Quenching", "distfunc": "Trapezoid"})
    param(quenching, "A", 100, 50, 100, False, "Quenching")
    param(quenching, "B", 50, 50, 50, False, "Quenching")
    competition = ET.SubElement(inp, "Competition")
    param(competition, "Window", 500, 500, 500, False, "Window")
    param(competition, "Shift", 50, 50, 50, False, "Promoter")
    promoters = ET.SubElement(inp, "Promoters")
    prom = ET.SubElement(promoters, "Promoter", {"name": "synCRE", "function": "Arrhenius2"})
    param(prom, "Q", 1, 1, 1, False)
    for name in ("Rmax", "Theta"):
        p = promoter[name]
        param(prom, name, p["value"], p["lim_low"], p["lim_high"], False if fixed else p["anneal"])

    tf_root = ET.SubElement(inp, "TFs")
    for tf_name, params in tfs.items():
        if tf_name not in pwms:
            raise ValueError(f"Missing PWM for {tf_name}")
        tf = ET.SubElement(tf_root, "TF", {"name": tf_name, "bsize": str(len(pwms[tf_name])), "include": "true"})
        for pname in ("kmax", "threshold", "lambda"):
            p = params[pname]
            param(tf, pname, p["value"], p["lim_low"], p["lim_high"], False if fixed else p["anneal"])
        coefs = ET.SubElement(tf, "Coefficients")
        p = params["coef"]
        param(coefs, "coef", p["value"], p["lim_low"], p["lim_high"], False if fixed else p["anneal"], "Promoter")
        pwm = ET.SubElement(tf, "PWM", {"type": "PSSM", "source": "Kang24_Table_S2"})
        for position in pwms[tf_name]:
            ET.SubElement(pwm, "position").text = "; ".join(fmt(v) for v in position)

    ET.SubElement(inp, "Interactions")
    genes = ET.SubElement(inp, "Genes")
    source = ET.SubElement(genes, "Source", {"name": "Kang24", "file": "local", "type": "local"})
    for _, row in seqs.iterrows():
        ET.SubElement(
            source,
            "Gene",
            {
                "name": row["xml_name"],
                "header": str(row["name"]),
                "promoter": "synCRE",
                "sequence": str(row["sequence"]).upper(),
                "include": "true",
            },
        )

    ET.SubElement(inp, "ScaleFactors")
    rate = ET.SubElement(inp, "RateData", {"row": "ID", "col": "gene"})
    rate_attrs = {"ID": "expr"}
    for _, row in seqs.iterrows():
        rate_attrs[row["xml_name"]] = fmt(row["expression level"])
    ET.SubElement(rate, "TableRow", rate_attrs)

    tfdata = ET.SubElement(inp, "TFData", {"row": "ID", "col": "TF"})
    tf_attrs = {"ID": "expr"}
    for tf_name in tfs:
        tf_attrs[tf_name] = "1"
    ET.SubElement(tfdata, "TableRow", tf_attrs)

    if mode := promoter.get("sheet"):
        root.insert(0, ET.Comment(f"Kang24 {mode}; model_num={fmt(promoter['model_num'])}; table_score_sse={fmt(promoter['score_sse'])}"))
    return ET.ElementTree(root)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", default="4TF_self", choices=["4TF_self", "7TF", "7TF_self"])
    parser.add_argument("--data-sheet", default="single-hit")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--model-column", type=int, help="0-based Excel column; default chooses lowest SSE column")
    parser.add_argument("--tf-names", help="Comma-separated subset of TFs to include, e.g. CREB1")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--seed", default="0")
    parser.add_argument("--init-loop", type=int)
    parser.add_argument("--trainable", action="store_true", help="Keep Table S1 tweak flags annealable; default fixes trained parameters")
    parser.add_argument("--competition", choices=["true", "false"], default="false", help="Use transcpp promoter/enhancer competition mode. Kang24 fits should use false.")
    parser.add_argument("--self-competition", choices=["true", "false"], default="true", help="Allow same-TF overlapping sites to self-compete.")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    promoter, tfs = read_model_parameters(SUPPLEMENTARY / "mmc2.xlsx", args.model, args.model_column)
    if args.tf_names:
        keep = [name.strip() for name in args.tf_names.split(",") if name.strip()]
        missing = [name for name in keep if name not in tfs]
        if missing:
            raise ValueError(f"Requested TFs are not in {args.model}: {missing}")
        tfs = {name: tfs[name] for name in keep}
    pwms = read_pwms(SUPPLEMENTARY / "mmc3.xlsx")
    seqs = read_sequences(SUPPLEMENTARY / "mmc4.xlsx", args.data_sheet, args.limit)
    tree = build_xml(
        promoter,
        tfs,
        pwms,
        seqs,
        smoke=args.smoke,
        fixed=not args.trainable,
        seed=args.seed,
        init_loop=args.init_loop,
        competition=args.competition == "true",
        self_competition=args.self_competition == "true",
    )
    ET.indent(tree, space="  ")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    tree.write(args.out, encoding="utf-8", xml_declaration=True)
    print(f"wrote {args.out}")
    print(f"model={args.model} model_num={fmt(promoter['model_num'])} table_score_sse={fmt(promoter['score_sse'])} genes={len(seqs)} tfs={len(tfs)}")


if __name__ == "__main__":
    main()
