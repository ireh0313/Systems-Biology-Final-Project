#!/usr/bin/env python3
"""Create Figure 5-style functional TFBS analysis from the fitted 4TF model."""

from __future__ import annotations

import argparse
import math
import re
import subprocess
from io import StringIO
from pathlib import Path
import xml.etree.ElementTree as ET

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TRANSC = ROOT / "work" / "repos" / "transcpp"
SUPPLEMENTARY = ROOT / "work" / "kang24" / "supplementary"
DEFAULT_XML = ROOT / "outputs" / "fullfit_final" / "kang24_4tf_full_seed0.xml"
DEFAULT_OUT = ROOT / "outputs" / "fig5_analysis"

TF_COLORS = {
    "ATF1": "#4c78a8",
    "ATF7": "#72b7b2",
    "CREB1": "#9a4f34",
    "CREM": "#b279a2",
}
BASE_COLORS = {"A": "#2f6f73", "T": "#9a4f34"}
CRE_REGIONS = {
    "CRE1": (31, 38),
    "CRE4": (89, 96),
}


def valid_xml_name(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]", "_", str(name))
    if not re.match(r"^[A-Za-z_]", cleaned):
        cleaned = f"g_{cleaned}"
    return cleaned


def read_observed() -> pd.DataFrame:
    df = pd.read_excel(SUPPLEMENTARY / "mmc4.xlsx", sheet_name="single-hit")
    df = df.dropna(subset=["name", "sequence", "expression level"]).copy()
    df["xml_name"] = [valid_xml_name(v) for v in df["name"]]
    extracted = df["name"].astype(str).str.extract(r"scanmut_single_pos_(\d+)_([ACGT])")
    df["position"] = pd.to_numeric(extracted[0], errors="coerce")
    df["mut_base"] = extracted[1]
    return df


def run_unfold(xml: Path, flag: str, gene: str | None = None) -> str:
    cmd = [str(TRANSC / "unfold"), "-i", str(xml.resolve()), "-s", "Output", flag]
    if gene:
        cmd.extend(["--gene", gene])
    result = subprocess.run(cmd, cwd=TRANSC, check=True, text=True, capture_output=True)
    return result.stdout


def read_rate(xml: Path) -> pd.DataFrame:
    text = run_unfold(xml, "--rate")
    row = pd.read_csv(StringIO(text), sep=r"\s+").iloc[0].drop(labels=["id"], errors="ignore")
    out = row.rename_axis("xml_name").reset_index(name="predicted")
    out["predicted"] = pd.to_numeric(out["predicted"])
    return out


def parse_unfold_table(text: str) -> pd.DataFrame:
    lines = [line for line in text.splitlines() if line.strip()]
    if len(lines) < 2:
        return pd.DataFrame()
    first_fields = lines[0].split()
    if first_fields and first_fields[0] == "id":
        return pd.read_csv(StringIO("\n".join(lines)), sep=r"\s+")
    # For --sites, first non-empty line is the gene name; the remaining text is a whitespace table.
    return pd.read_csv(StringIO("\n".join(lines[1:])), sep=r"\s+")


def sites_for_gene(xml: Path, gene: str) -> pd.DataFrame:
    sites = parse_unfold_table(run_unfold(xml, "--sites", gene))
    occ = parse_unfold_table(run_unfold(xml, "--occupancy", gene))
    if sites.empty or occ.empty:
        return pd.DataFrame()
    occ_long = occ.drop(columns=["id"]).T.reset_index()
    occ_long.columns = ["site_id", "occupancy"]
    tf_names = sorted(sites["name"].unique(), key=len, reverse=True)

    def split_site_id(site_id: str) -> tuple[str | None, int | None]:
        for tf_name in tf_names:
            if site_id.startswith(tf_name):
                suffix = site_id[len(tf_name) :]
                if suffix.isdigit():
                    return tf_name, int(suffix)
        return None, None

    split = occ_long["site_id"].apply(split_site_id)
    occ_long["tf"] = [item[0] for item in split]
    occ_long["index"] = [item[1] for item in split]
    merged = sites.merge(occ_long, left_on=["name", "index"], right_on=["tf", "index"], how="left")
    merged["seq_start"] = merged["start"] + 194
    merged["seq_end"] = merged["end"] + 194
    merged["center"] = (merged["seq_start"] + merged["seq_end"]) / 2
    return merged


def coefficients(xml: Path) -> dict[str, float]:
    root = ET.parse(xml).getroot()
    out: dict[str, float] = {}
    section = root.find("Output")
    if section is None:
        section = root.find("Input")
    if section is None:
        return out
    for tf in section.findall("./TFs/TF"):
        name = tf.attrib["name"]
        coef = tf.find("./Coefficients/coef")
        if coef is not None:
            out[name] = float(coef.attrib["value"])
    return out


def add_delta(df: pd.DataFrame) -> pd.DataFrame:
    wt = df.loc[df["name"].eq("synCRE_Promega_0")].iloc[0]
    out = df.copy()
    out["observed_delta"] = (out["expression level"] / wt["expression level"]).apply(lambda x: math.log2(x) if x > 0 else math.nan)
    out["predicted_delta"] = (out["predicted"] / wt["predicted"]).apply(lambda x: math.log2(x) if x > 0 else math.nan)
    return out


def svg_escape(text: object) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def scale(v: float, src_min: float, src_max: float, dst_min: float, dst_max: float) -> float:
    if src_max == src_min:
        return (dst_min + dst_max) / 2
    return dst_min + (v - src_min) / (src_max - src_min) * (dst_max - dst_min)


def write_svg(path: Path, width: int, height: int, body: list[str]) -> None:
    path.write_text(
        "\n".join([
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
            '<rect width="100%" height="100%" fill="#ffffff"/>',
            *body,
            "</svg>",
        ]),
        encoding="utf-8",
    )


def line_points(df: pd.DataFrame, xcol: str, ycol: str, xlo: float, xhi: float, ylo: float, yhi: float, left: int, top: int, plot_w: int, plot_h: int) -> str:
    pairs = []
    for row in df.sort_values(xcol).itertuples(index=False):
        x = scale(float(getattr(row, xcol)), xlo, xhi, left, left + plot_w)
        y = scale(float(getattr(row, ycol)), ylo, yhi, top + plot_h, top)
        pairs.append(f"{x:.2f},{y:.2f}")
    return " ".join(pairs)


def make_at_profile(df: pd.DataFrame, out: Path) -> None:
    plot = df[df["mut_base"].isin(["A", "T"])].dropna(subset=["position"]).copy()
    width, height = 1080, 620
    left, right, top, bottom = 92, 45, 72, 82
    plot_w, plot_h = width - left - right, height - top - bottom
    xlo, xhi = 0, 86
    ylo = min(plot["observed_delta"].min(), plot["predicted_delta"].min())
    yhi = max(plot["observed_delta"].max(), plot["predicted_delta"].max())
    ypad = max(0.1, (yhi - ylo) * 0.08)
    ylo, yhi = ylo - ypad, yhi + ypad
    body = [
        f'<text x="{width/2}" y="34" text-anchor="middle" font-family="Arial" font-size="22" font-weight="700">Figure 5A/G: A/T substitution activity profiles</text>',
        f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="#fbfbfb" stroke="#222"/>',
    ]
    for pos, label in [(11, "CRE1"), (69, "CRE4")]:
        x = scale(pos, xlo, xhi, left, left + plot_w)
        body.append(f'<rect x="{x:.2f}" y="{top}" width="{scale(pos+7, xlo, xhi, left, left+plot_w)-x:.2f}" height="{plot_h}" fill="#f2e7d8" opacity="0.65"/>')
        body.append(f'<text x="{x+34:.2f}" y="{top+18}" text-anchor="middle" font-family="Arial" font-size="12" fill="#6b4a29">{label}</text>')
    for t in [-2, -1, 0, 1]:
        y = scale(t, ylo, yhi, top + plot_h, top)
        body.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left+plot_w}" y2="{y:.2f}" stroke="#e1e1e1"/>')
        body.append(f'<text x="{left-12}" y="{y+4:.2f}" text-anchor="end" font-family="Arial" font-size="13">{t}</text>')
    for t in range(0, 87, 10):
        x = scale(t, xlo, xhi, left, left + plot_w)
        body.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top+plot_h}" stroke="#ececec"/>')
        body.append(f'<text x="{x:.2f}" y="{top+plot_h+25}" text-anchor="middle" font-family="Arial" font-size="13">{t}</text>')
    dash = {"A": "", "T": "stroke-dasharray=\"6 5\""}
    for base in ["A", "T"]:
        sub = plot[plot["mut_base"].eq(base)]
        color = BASE_COLORS[base]
        obs = line_points(sub, "position", "observed_delta", xlo, xhi, ylo, yhi, left, top, plot_w, plot_h)
        pred = line_points(sub, "position", "predicted_delta", xlo, xhi, ylo, yhi, left, top, plot_w, plot_h)
        body.append(f'<polyline fill="none" stroke="{color}" stroke-width="2.5" opacity="0.55" points="{obs}"/>')
        body.append(f'<polyline fill="none" stroke="{color}" stroke-width="3.5" {dash[base]} points="{pred}"/>')
    body.extend([
        f'<text x="{left+plot_w/2}" y="{height-26}" text-anchor="middle" font-family="Arial" font-size="16">Single-hit mutation position</text>',
        f'<text x="25" y="{top+plot_h/2}" text-anchor="middle" font-family="Arial" font-size="16" transform="rotate(-90 25 {top+plot_h/2})">Delta activity, log2(variant / WT)</text>',
        f'<text x="{width-288}" y="42" font-family="Arial" font-size="13" fill="{BASE_COLORS["A"]}">A substitution</text>',
        f'<text x="{width-166}" y="42" font-family="Arial" font-size="13" fill="{BASE_COLORS["T"]}">T substitution</text>',
        f'<text x="{width-288}" y="60" font-family="Arial" font-size="12" fill="#555">thin: observed, thick/dashed: model</text>',
    ])
    write_svg(out, width, height, body)


def region_label(seq_start: float, seq_end: float) -> str | None:
    for name, (lo, hi) in CRE_REGIONS.items():
        if seq_start <= hi and seq_end >= lo:
            return name
    return None


def contribution_table(xml: Path, genes: list[str], coefs: dict[str, float]) -> pd.DataFrame:
    rows = []
    for gene in genes:
        sites = sites_for_gene(xml, gene)
        if sites.empty:
            continue
        for row in sites.itertuples(index=False):
            region = region_label(float(row.seq_start), float(row.seq_end))
            if not region:
                continue
            coef = coefs.get(row.name, math.nan)
            rows.append({
                "gene": gene,
                "tf": row.name,
                "region": region,
                "site_index": row.index,
                "seq_start": row.seq_start,
                "seq_end": row.seq_end,
                "score": row.score,
                "occupancy": row.occupancy,
                "coef": coef,
                "activation_contribution": row.occupancy * coef,
            })
    return pd.DataFrame(rows)


def make_activation_coefficients(coefs: dict[str, float], out: Path) -> None:
    width, height = 620, 390
    left, top, plot_w, plot_h = 88, 60, 470, 250
    maxv = max(coefs.values())
    body = [
        f'<text x="{width/2}" y="32" text-anchor="middle" font-family="Arial" font-size="21" font-weight="700">Figure 5E: Activation coefficients</text>',
    ]
    for i, (tf, value) in enumerate(sorted(coefs.items())):
        x = left + i * (plot_w / len(coefs)) + 18
        w = plot_w / len(coefs) - 34
        h = scale(value, 0, maxv, 0, plot_h)
        y = top + plot_h - h
        body.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" fill="{TF_COLORS.get(tf, "#777")}"/>')
        body.append(f'<text x="{x+w/2:.2f}" y="{top+plot_h+24}" text-anchor="middle" font-family="Arial" font-size="13">{tf}</text>')
        body.append(f'<text x="{x+w/2:.2f}" y="{y-7:.2f}" text-anchor="middle" font-family="Arial" font-size="12">{value:.2g}</text>')
    body.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="#222"/>')
    body.append(f'<text x="24" y="{top+plot_h/2}" text-anchor="middle" font-family="Arial" font-size="15" transform="rotate(-90 24 {top+plot_h/2})">Coefficient</text>')
    write_svg(out, width, height, body)


def make_contribution_bars(contrib: pd.DataFrame, out: Path) -> None:
    selected = ["synCRE_Promega_0", "scanmut_single_pos_14_A", "scanmut_single_pos_72_A", "scanmut_single_pos_14_T", "scanmut_single_pos_72_T"]
    summary = contrib.groupby(["gene", "region", "tf"], as_index=False)["activation_contribution"].sum()
    summary = summary[summary["gene"].isin(selected)].copy()
    width, height = 1120, 560
    left, top, plot_w, plot_h = 95, 70, 970, 370
    groups = [
        ("CRE1 WT", "synCRE_Promega_0", "CRE1"),
        ("CRE1 C14A", "scanmut_single_pos_14_A", "CRE1"),
        ("CRE4 WT", "synCRE_Promega_0", "CRE4"),
        ("CRE4 C72A", "scanmut_single_pos_72_A", "CRE4"),
        ("CRE1 C14T", "scanmut_single_pos_14_T", "CRE1"),
        ("CRE4 C72T", "scanmut_single_pos_72_T", "CRE4"),
    ]
    maxv = summary.groupby(["gene", "region"])["activation_contribution"].sum().max()
    bar_w = plot_w / len(groups) * 0.58
    body = [
        f'<text x="{width/2}" y="34" text-anchor="middle" font-family="Arial" font-size="22" font-weight="700">Figure 5F/K: CRE1/CRE4 activation contribution</text>',
        f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="#222"/>',
    ]
    for t in [0, maxv / 2, maxv]:
        y = scale(t, 0, maxv, top + plot_h, top)
        body.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left+plot_w}" y2="{y:.2f}" stroke="#e5e5e5"/>')
        body.append(f'<text x="{left-10}" y="{y+4:.2f}" text-anchor="end" font-family="Arial" font-size="12">{t:.2f}</text>')
    for i, (label, gene, region) in enumerate(groups):
        x = left + i * (plot_w / len(groups)) + 25
        bottom = top + plot_h
        rows = summary[(summary["gene"].eq(gene)) & (summary["region"].eq(region))].sort_values("tf")
        for row in rows.itertuples(index=False):
            h = scale(float(row.activation_contribution), 0, maxv, 0, plot_h)
            bottom -= h
            body.append(f'<rect x="{x:.2f}" y="{bottom:.2f}" width="{bar_w:.2f}" height="{h:.2f}" fill="{TF_COLORS.get(row.tf, "#777")}" stroke="#fff" stroke-width="0.5"/>')
        body.append(f'<text x="{x+bar_w/2:.2f}" y="{top+plot_h+24}" text-anchor="middle" font-family="Arial" font-size="12">{svg_escape(label)}</text>')
    lx = width - 205
    for i, tf in enumerate(sorted(TF_COLORS)):
        body.append(f'<rect x="{lx}" y="{82+i*22}" width="13" height="13" fill="{TF_COLORS[tf]}"/>')
        body.append(f'<text x="{lx+20}" y="{93+i*22}" font-family="Arial" font-size="13">{tf}</text>')
    body.append(f'<text x="26" y="{top+plot_h/2}" text-anchor="middle" font-family="Arial" font-size="15" transform="rotate(-90 26 {top+plot_h/2})">sum(coef x occupancy)</text>')
    write_svg(out, width, height, body)


def make_site_maps(xml: Path, genes: list[str], coefs: dict[str, float], out: Path) -> pd.DataFrame:
    all_sites = []
    for gene in genes:
        sites = sites_for_gene(xml, gene)
        sites["gene"] = gene
        sites["coef"] = sites["name"].map(coefs)
        sites["activation_contribution"] = sites["occupancy"] * sites["coef"]
        all_sites.append(sites)
    sites_df = pd.concat(all_sites, ignore_index=True)
    focus = sites_df[(sites_df["seq_start"] <= 100) & (sites_df["seq_end"] >= 25)].copy()
    width = 1160
    row_h = 118
    height = 70 + row_h * len(genes) + 36
    left, right = 88, 35
    plot_w = width - left - right
    xlo, xhi = 25, 100
    body = [
        f'<text x="{width/2}" y="34" text-anchor="middle" font-family="Arial" font-size="22" font-weight="700">Figure 5B/C/D/I/J: TFBS occupancy around CRE1 and CRE4</text>',
    ]
    for gi, gene in enumerate(genes):
        y0 = 68 + gi * row_h
        label = gene.replace("scanmut_single_pos_", "pos ")
        body.append(f'<text x="{left-12}" y="{y0+42}" text-anchor="end" font-family="Arial" font-size="13">{svg_escape(label)}</text>')
        body.append(f'<line x1="{left}" y1="{y0+42}" x2="{left+plot_w}" y2="{y0+42}" stroke="#ddd"/>')
        for region, (lo, hi) in CRE_REGIONS.items():
            x = scale(lo, xlo, xhi, left, left + plot_w)
            w = scale(hi + 1, xlo, xhi, left, left + plot_w) - x
            body.append(f'<rect x="{x:.2f}" y="{y0+4}" width="{w:.2f}" height="78" fill="#f2e7d8" opacity="0.65"/>')
            body.append(f'<text x="{x+w/2:.2f}" y="{y0+16}" text-anchor="middle" font-family="Arial" font-size="11" fill="#6b4a29">{region}</text>')
        sub = focus[focus["gene"].eq(gene)].sort_values(["name", "seq_start"])
        lanes = {tf: i for i, tf in enumerate(sorted(TF_COLORS))}
        for row in sub.itertuples(index=False):
            lane = lanes.get(row.name, 0)
            x = scale(float(row.seq_start), xlo, xhi, left, left + plot_w)
            w = max(3, scale(float(row.seq_end), xlo, xhi, left, left + plot_w) - x)
            y = y0 + 23 + lane * 15
            opacity = max(0.15, min(1.0, float(row.occupancy) / 0.30))
            body.append(f'<rect x="{x:.2f}" y="{y}" width="{w:.2f}" height="10" fill="{TF_COLORS.get(row.name, "#777")}" opacity="{opacity:.2f}"/>')
        for pos in [31, 38, 89, 96]:
            x = scale(pos, xlo, xhi, left, left + plot_w)
            body.append(f'<line x1="{x:.2f}" y1="{y0+4}" x2="{x:.2f}" y2="{y0+86}" stroke="#c7a66f" stroke-dasharray="3 3"/>')
    body.append(f'<text x="{left+plot_w/2}" y="{height-10}" text-anchor="middle" font-family="Arial" font-size="15">Sequence coordinate in full construct; CRE1 = 31-38, CRE4 = 89-96</text>')
    write_svg(out, width, height, body)
    return sites_df


def write_summary(out: Path, corr: float, coefs: dict[str, float], key_rows: pd.DataFrame) -> None:
    lines = [
        "# Kang24 Figure 5-Style Functional Binding Site Analysis",
        "",
        "This analysis uses the completed 4TF full fitted model and `unfold` outputs for sites and fractional occupancy.",
        "",
        f"- 4TF model A/T substitution Pearson r: `{corr:.6f}`",
        "- Contribution proxy used here: `activation coefficient x fractional occupancy`.",
        "- CRE1 and CRE4 coordinates are mapped from the single-hit scan: CRE1 `31-38`, CRE4 `89-96` in the full construct, corresponding to scan positions `11-18` and `69-76`.",
        "",
        "## Activation coefficients",
        "",
    ]
    for tf, coef in sorted(coefs.items()):
        lines.append(f"- {tf}: `{coef:.6g}`")
    lines.extend(["", "## Key variants", "", "| Variant | Observed delta | Predicted delta |", "|---|---:|---:|"])
    for row in key_rows.itertuples(index=False):
        lines.append(f"| {row.name} | `{row.observed_delta:.6f}` | `{row.predicted_delta:.6f}` |")
    lines.extend([
        "",
        "## Outputs",
        "",
        "- `kang24_fig5_AT_substitution_profile.svg`",
        "- `kang24_fig5_TFBS_occupancy_maps.svg`",
        "- `kang24_fig5_activation_coefficients.svg`",
        "- `kang24_fig5_CRE1_CRE4_contribution_bars.svg`",
        "- `kang24_fig5_sites_occupancy_contribution.csv`",
        "- `kang24_fig5_AT_substitution_table.csv`",
    ])
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--xml", type=Path, default=DEFAULT_XML)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    observed = read_observed()
    merged = add_delta(observed.merge(read_rate(args.xml), on="xml_name", how="inner"))
    at = merged[merged["mut_base"].isin(["A", "T"])].dropna(subset=["observed_delta", "predicted_delta"]).copy()
    at.to_csv(args.outdir / "kang24_fig5_AT_substitution_table.csv", index=False)
    make_at_profile(merged, args.outdir / "kang24_fig5_AT_substitution_profile.svg")

    coefs = coefficients(args.xml)
    genes = [
        "synCRE_Promega_0",
        "scanmut_single_pos_14_A",
        "scanmut_single_pos_72_A",
        "scanmut_single_pos_14_T",
        "scanmut_single_pos_72_T",
    ]
    sites_df = make_site_maps(args.xml, genes, coefs, args.outdir / "kang24_fig5_TFBS_occupancy_maps.svg")
    sites_df.to_csv(args.outdir / "kang24_fig5_sites_occupancy_contribution.csv", index=False)
    contrib = contribution_table(args.xml, genes, coefs)
    contrib.to_csv(args.outdir / "kang24_fig5_CRE_region_contribution_table.csv", index=False)
    make_activation_coefficients(coefs, args.outdir / "kang24_fig5_activation_coefficients.svg")
    make_contribution_bars(contrib, args.outdir / "kang24_fig5_CRE1_CRE4_contribution_bars.svg")
    key_rows = merged[merged["name"].isin(genes)].copy()
    write_summary(args.outdir / "kang24_fig5_summary.md", float(at["observed_delta"].corr(at["predicted_delta"])), coefs, key_rows)
    print(f"wrote {args.outdir}")
    print(f"AT Pearson r: {at['observed_delta'].corr(at['predicted_delta']):.6f}")


if __name__ == "__main__":
    main()
