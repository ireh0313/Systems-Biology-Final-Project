#!/usr/bin/env python3
"""Create Kang24 Figure 4J-style multi-hit prediction from the fitted 4TF model."""

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
DEFAULT_MODEL = ROOT / "outputs" / "fullfit_final" / "kang24_4tf_full_seed0.xml"
DEFAULT_OUT = ROOT / "outputs" / "fig4j_multihit"


def valid_xml_name(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]", "_", str(name))
    if not re.match(r"^[A-Za-z_]", cleaned):
        cleaned = f"g_{cleaned}"
    return cleaned


def fmt(value: object) -> str:
    if isinstance(value, str):
        return value
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return f"{float(value):.12g}"


def read_multihit(sheet: str) -> pd.DataFrame:
    df = pd.read_excel(SUPPLEMENTARY / "mmc4.xlsx", sheet_name=sheet)
    df = df.dropna(subset=["name", "sequence", "expression level"]).copy()
    df["xml_name"] = [valid_xml_name(v) for v in df["name"]]
    if df["xml_name"].duplicated().any():
        counts: dict[str, int] = {}
        names = []
        for name in df["xml_name"]:
            counts[name] = counts.get(name, 0) + 1
            names.append(name if counts[name] == 1 else f"{name}_{counts[name]}")
        df["xml_name"] = names
    return df


def replace_genes_and_rate(section: ET.Element, seqs: pd.DataFrame) -> None:
    genes = section.find("Genes")
    if genes is None:
        genes = ET.SubElement(section, "Genes")
    for child in list(genes):
        genes.remove(child)
    source = ET.SubElement(genes, "Source", {"name": "Kang24_multi_hit", "file": "local", "type": "local"})
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

    ratedata = section.find("RateData")
    if ratedata is None:
        ratedata = ET.SubElement(section, "RateData", {"row": "ID", "col": "gene"})
    else:
        ratedata.attrib.update({"row": "ID", "col": "gene"})
    for child in list(ratedata):
        ratedata.remove(child)
    attrs = {"ID": "expr"}
    for _, row in seqs.iterrows():
        attrs[row["xml_name"]] = fmt(row["expression level"])
    ET.SubElement(ratedata, "TableRow", attrs)


def build_multihit_xml(model_xml: Path, seqs: pd.DataFrame, out_xml: Path) -> None:
    tree = ET.parse(model_xml)
    root = tree.getroot()
    for section_name in ["Input", "Output"]:
        section = root.find(section_name)
        if section is not None:
            replace_genes_and_rate(section, seqs)
    root.insert(0, ET.Comment("Figure 4J multi-hit prediction XML generated from fitted 4TF full model; parameters unchanged"))
    ET.indent(tree, space="  ")
    out_xml.parent.mkdir(parents=True, exist_ok=True)
    tree.write(out_xml, encoding="utf-8", xml_declaration=True)


def run_rate(xml: Path) -> pd.DataFrame:
    cmd = [str(TRANSC / "unfold"), "-i", str(xml.resolve()), "-s", "Output", "--rate"]
    result = subprocess.run(cmd, cwd=TRANSC, check=True, text=True, capture_output=True)
    row = pd.read_csv(StringIO(result.stdout), sep=r"\s+").iloc[0].drop(labels=["id"], errors="ignore")
    out = row.rename_axis("xml_name").reset_index(name="predicted")
    out["predicted"] = pd.to_numeric(out["predicted"])
    return out


def add_delta(df: pd.DataFrame) -> pd.DataFrame:
    wt = df.loc[df["name"].astype(str).str.contains("synCRE_Promega_0", regex=False)]
    if wt.empty:
        raise ValueError("Could not find multi-hit WT row")
    observed_wt = float(wt["expression level"].iloc[0])
    predicted_wt = float(wt["predicted"].iloc[0])
    out = df.copy()
    out["observed_delta"] = (out["expression level"] / observed_wt).apply(lambda x: math.log2(x) if x > 0 else math.nan)
    out["predicted_delta"] = (out["predicted"] / predicted_wt).apply(lambda x: math.log2(x) if x > 0 else math.nan)
    return out.dropna(subset=["observed_delta", "predicted_delta"])


def scale(value: float, src_min: float, src_max: float, dst_min: float, dst_max: float) -> float:
    if src_max == src_min:
        return (dst_min + dst_max) / 2
    return dst_min + (value - src_min) / (src_max - src_min) * (dst_max - dst_min)


def ticks(lo: float, hi: float, n: int = 6) -> list[float]:
    if lo == hi:
        return [lo]
    step = (hi - lo) / (n - 1)
    return [lo + i * step for i in range(n)]


def write_svg(path: Path, width: int, height: int, body: list[str]) -> None:
    path.write_text(
        "\n".join(
            [
                f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
                '<rect width="100%" height="100%" fill="#ffffff"/>',
                *body,
                "</svg>",
            ]
        ),
        encoding="utf-8",
    )


def make_scatter(df: pd.DataFrame, out: Path) -> tuple[float, float]:
    corr = float(df["observed_delta"].corr(df["predicted_delta"]))
    rmse = math.sqrt(float(((df["observed_delta"] - df["predicted_delta"]) ** 2).mean()))
    width, height = 880, 720
    left, right, top, bottom = 105, 42, 92, 92
    plot_w, plot_h = width - left - right, height - top - bottom
    lo = min(float(df["observed_delta"].min()), float(df["predicted_delta"].min()))
    hi = max(float(df["observed_delta"].max()), float(df["predicted_delta"].max()))
    pad = max(0.1, (hi - lo) * 0.06)
    lo, hi = lo - pad, hi + pad
    body = [
        f'<text x="{width/2}" y="34" text-anchor="middle" font-family="Arial" font-size="22" font-weight="700">Figure 4J: Multi-hit prediction with 4TF model</text>',
        f'<text x="{width/2}" y="62" text-anchor="middle" font-family="Arial" font-size="15" fill="#555">Pearson r = {corr:.3f}, RMSE = {rmse:.3f}, n = {len(df)}</text>',
        f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="#fbfbfb" stroke="#222"/>',
    ]
    for t in ticks(lo, hi):
        x = scale(t, lo, hi, left, left + plot_w)
        y = scale(t, lo, hi, top + plot_h, top)
        body.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top+plot_h}" stroke="#e1e1e1"/>')
        body.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left+plot_w}" y2="{y:.2f}" stroke="#e1e1e1"/>')
        body.append(f'<text x="{x:.2f}" y="{top+plot_h+25}" text-anchor="middle" font-family="Arial" font-size="13">{t:.1f}</text>')
        body.append(f'<text x="{left-12}" y="{y+4:.2f}" text-anchor="end" font-family="Arial" font-size="13">{t:.1f}</text>')
    body.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top}" stroke="#9a4f34" stroke-width="2" stroke-dasharray="7 6"/>')
    for row in df.itertuples(index=False):
        x = scale(float(row.observed_delta), lo, hi, left, left + plot_w)
        y = scale(float(row.predicted_delta), lo, hi, top + plot_h, top)
        body.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3.2" fill="#2f6f73" opacity="0.58"/>')
    body.extend(
        [
            f'<text x="{left+plot_w/2}" y="{height-28}" text-anchor="middle" font-family="Arial" font-size="16">Observed multi-hit Δactivity, log2(variant / WT)</text>',
            f'<text x="25" y="{top+plot_h/2}" text-anchor="middle" font-family="Arial" font-size="16" transform="rotate(-90 25 {top+plot_h/2})">Predicted multi-hit Δactivity, log2(variant / WT)</text>',
        ]
    )
    write_svg(out, width, height, body)
    return corr, rmse


def write_summary(path: Path, corr: float, rmse: float, n: int, xml: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "# Kang24 Figure 4J Multi-Hit Prediction",
                "",
                "This reproduces the Figure 4J-style multi-hit validation using the completed 4TF full fitted model.",
                "",
                f"- Model XML with multi-hit genes: `{xml.name}`",
                f"- Multi-hit variants included: `{n}`",
                f"- Pearson r: `{corr:.6f}`",
                f"- RMSE: `{rmse:.6f}`",
                "",
                "The fitted 4TF parameters were not retrained on multi-hit data. Only the gene sequence table and observed RateData were replaced with the `multi-hit` sheet from Kang24 Table S3.",
                "",
                "## Outputs",
                "",
                "- `kang24_fig4j_4tf_multihit_scatter.svg`",
                "- `kang24_fig4j_4tf_multihit_observed_predicted.csv`",
                "- `kang24_4tf_full_seed0_multihit.xml`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-xml", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--sheet", default="multi-hit")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    seqs = read_multihit(args.sheet)
    multihit_xml = args.outdir / "kang24_4tf_full_seed0_multihit.xml"
    build_multihit_xml(args.model_xml, seqs, multihit_xml)
    predicted = run_rate(multihit_xml)
    merged = add_delta(seqs.merge(predicted, on="xml_name", how="inner"))
    merged.to_csv(args.outdir / "kang24_fig4j_4tf_multihit_observed_predicted.csv", index=False)
    corr, rmse = make_scatter(merged, args.outdir / "kang24_fig4j_4tf_multihit_scatter.svg")
    write_summary(args.outdir / "kang24_fig4j_summary.md", corr, rmse, len(merged), multihit_xml)
    print(f"wrote {args.outdir}")
    print(f"Pearson r: {corr:.6f}")
    print(f"RMSE: {rmse:.6f}")
    print(f"n: {len(merged)}")


if __name__ == "__main__":
    main()
