#!/usr/bin/env python3
"""Create Kang24 observed-vs-predicted figures from a transcpp XML file."""

from __future__ import annotations

import argparse
import re
import subprocess
from io import StringIO
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TRANSC = ROOT / "work" / "repos" / "transcpp"
SUPPLEMENTARY = ROOT / "work" / "kang24" / "supplementary"


def valid_xml_name(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]", "_", str(name))
    if not re.match(r"^[A-Za-z_]", cleaned):
        cleaned = f"g_{cleaned}"
    return cleaned


def run_unfold_rate(xml: Path, section: str) -> pd.DataFrame:
    cmd = [str(TRANSC / "unfold"), "-i", str(xml.resolve()), "-s", section, "--rate"]
    result = subprocess.run(cmd, cwd=TRANSC, check=True, text=True, capture_output=True)
    return pd.read_csv(StringIO(result.stdout), sep=r"\s+")


def read_observed(sheet: str) -> pd.DataFrame:
    df = pd.read_excel(SUPPLEMENTARY / "mmc4.xlsx", sheet_name=sheet)
    df = df.dropna(subset=["name", "sequence", "expression level"]).copy()
    df["xml_name"] = [valid_xml_name(v) for v in df["name"]]
    if df["xml_name"].duplicated().any():
        counts = {}
        names = []
        for name in df["xml_name"]:
            counts[name] = counts.get(name, 0) + 1
            names.append(name if counts[name] == 1 else f"{name}_{counts[name]}")
        df["xml_name"] = names
    return df


def predictions_to_long(rate_df: pd.DataFrame) -> pd.DataFrame:
    row = rate_df.iloc[0].drop(labels=["id"], errors="ignore")
    out = row.rename_axis("xml_name").reset_index(name="predicted")
    out["predicted"] = pd.to_numeric(out["predicted"])
    return out


def add_scan_position(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    extracted = out["name"].astype(str).str.extract(r"scanmut_single_pos_(\d+)_([ACGT])")
    out["position"] = pd.to_numeric(extracted[0], errors="coerce")
    out["mut_base"] = extracted[1]
    return out


def svg_escape(text: object) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def scale(value: float, src_min: float, src_max: float, dst_min: float, dst_max: float) -> float:
    if src_max == src_min:
        return (dst_min + dst_max) / 2
    return dst_min + (value - src_min) * (dst_max - dst_min) / (src_max - src_min)


def ticks(lo: float, hi: float, n: int = 5) -> list[float]:
    if hi == lo:
        return [lo]
    step = (hi - lo) / (n - 1)
    return [lo + i * step for i in range(n)]


def write_svg(out: Path, width: int, height: int, body: list[str]) -> None:
    out.write_text(
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


def add_delta_activity(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    wt = out.loc[out["name"].astype(str).str.contains("synCRE_Promega_0", regex=False)]
    if wt.empty:
        observed_wt = float(out["expression level"].iloc[0])
        predicted_wt = float(out["predicted"].iloc[0])
    else:
        observed_wt = float(wt["expression level"].iloc[0])
        predicted_wt = float(wt["predicted"].iloc[0])
    out["observed_delta"] = (out["expression level"] / observed_wt).apply(lambda x: pd.NA if x <= 0 else __import__("math").log2(x))
    out["predicted_delta"] = (out["predicted"] / predicted_wt).apply(lambda x: pd.NA if x <= 0 else __import__("math").log2(x))
    return out.dropna(subset=["observed_delta", "predicted_delta"])


def make_scatter(df: pd.DataFrame, out: Path, title: str, xcol: str, ycol: str, xlabel: str, ylabel: str) -> None:
    corr = df[xcol].corr(df[ycol], method="pearson")
    width, height = 900, 720
    left, right, top, bottom = 105, 40, 92, 95
    plot_w, plot_h = width - left - right, height - top - bottom
    lo = min(df[xcol].min(), df[ycol].min())
    hi = max(df[xcol].max(), df[ycol].max())
    pad = (hi - lo) * 0.06 if hi > lo else 1
    lo, hi = lo - pad, hi + pad
    body = [
        f'<text x="{width/2}" y="34" text-anchor="middle" font-family="Arial" font-size="22" font-weight="700">{svg_escape(title)}</text>',
        f'<text x="{width/2}" y="62" text-anchor="middle" font-family="Arial" font-size="16" fill="#555">Pearson r = {corr:.3f}, n = {len(df)}</text>',
        f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="#fbfbfb" stroke="#222" stroke-width="1"/>',
    ]
    for t in ticks(lo, hi):
        x = scale(t, lo, hi, left, left + plot_w)
        y = scale(t, lo, hi, top + plot_h, top)
        body.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top+plot_h}" stroke="#dddddd" stroke-width="1"/>')
        body.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left+plot_w}" y2="{y:.2f}" stroke="#dddddd" stroke-width="1"/>')
        body.append(f'<text x="{x:.2f}" y="{top+plot_h+26}" text-anchor="middle" font-family="Arial" font-size="13" fill="#444">{t:.0f}</text>')
        body.append(f'<text x="{left-13}" y="{y+4:.2f}" text-anchor="end" font-family="Arial" font-size="13" fill="#444">{t:.0f}</text>')
    body.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top}" stroke="#9a4f34" stroke-width="2" stroke-dasharray="7 6"/>')
    for _, row in df.iterrows():
        x = scale(float(row[xcol]), lo, hi, left, left + plot_w)
        y = scale(float(row[ycol]), lo, hi, top + plot_h, top)
        body.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="4.0" fill="#2f6f73" opacity="0.72"/>')
    body.extend(
        [
            f'<text x="{left+plot_w/2}" y="{height-32}" text-anchor="middle" font-family="Arial" font-size="16">{svg_escape(xlabel)}</text>',
            f'<text x="24" y="{top+plot_h/2}" text-anchor="middle" font-family="Arial" font-size="16" transform="rotate(-90 24 {top+plot_h/2})">{svg_escape(ylabel)}</text>',
        ]
    )
    write_svg(out, width, height, body)


def make_profile(df: pd.DataFrame, out: Path, title: str, observed_col: str, predicted_col: str, ylabel: str) -> None:
    scan = add_scan_position(df).dropna(subset=["position"]).copy()
    profile = (
        scan.groupby("position", as_index=False)
        .agg(observed=(observed_col, "mean"), predicted=(predicted_col, "mean"))
        .sort_values("position")
    )
    width, height = 1040, 560
    left, right, top, bottom = 95, 42, 72, 82
    plot_w, plot_h = width - left - right, height - top - bottom
    xlo, xhi = float(profile["position"].min()), float(profile["position"].max())
    ylo = min(float(profile["observed"].min()), float(profile["predicted"].min()))
    yhi = max(float(profile["observed"].max()), float(profile["predicted"].max()))
    ypad = (yhi - ylo) * 0.08 if yhi > ylo else 1
    ylo, yhi = ylo - ypad, yhi + ypad

    def points(col: str) -> str:
        pairs = []
        for _, row in profile.iterrows():
            x = scale(float(row["position"]), xlo, xhi, left, left + plot_w)
            y = scale(float(row[col]), ylo, yhi, top + plot_h, top)
            pairs.append(f"{x:.2f},{y:.2f}")
        return " ".join(pairs)

    body = [
        f'<text x="{width/2}" y="36" text-anchor="middle" font-family="Arial" font-size="22" font-weight="700">{svg_escape(title)}</text>',
        f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="#fbfbfb" stroke="#222" stroke-width="1"/>',
    ]
    for t in ticks(xlo, xhi, 8):
        x = scale(t, xlo, xhi, left, left + plot_w)
        body.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top+plot_h}" stroke="#e1e1e1" stroke-width="1"/>')
        body.append(f'<text x="{x:.2f}" y="{top+plot_h+25}" text-anchor="middle" font-family="Arial" font-size="13" fill="#444">{t:.0f}</text>')
    for t in ticks(ylo, yhi):
        y = scale(t, ylo, yhi, top + plot_h, top)
        body.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left+plot_w}" y2="{y:.2f}" stroke="#e1e1e1" stroke-width="1"/>')
        body.append(f'<text x="{left-12}" y="{y+4:.2f}" text-anchor="end" font-family="Arial" font-size="13" fill="#444">{t:.0f}</text>')
    body.append(f'<polyline fill="none" stroke="#2f6f73" stroke-width="3" points="{points("observed")}"/>')
    body.append(f'<polyline fill="none" stroke="#9a4f34" stroke-width="3" points="{points("predicted")}"/>')
    body.extend(
        [
            f'<line x1="{width-220}" y1="34" x2="{width-180}" y2="34" stroke="#2f6f73" stroke-width="3"/><text x="{width-170}" y="39" font-family="Arial" font-size="14">Observed</text>',
            f'<line x1="{width-220}" y1="58" x2="{width-180}" y2="58" stroke="#9a4f34" stroke-width="3"/><text x="{width-170}" y="63" font-family="Arial" font-size="14">Predicted</text>',
            f'<text x="{left+plot_w/2}" y="{height-26}" text-anchor="middle" font-family="Arial" font-size="16">Single-hit mutation position</text>',
            f'<text x="25" y="{top+plot_h/2}" text-anchor="middle" font-family="Arial" font-size="16" transform="rotate(-90 25 {top+plot_h/2})">{svg_escape(ylabel)}</text>',
        ]
    )
    write_svg(out, width, height, body)
    profile.to_csv(out.with_suffix(".csv"), index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--xml", type=Path, required=True)
    parser.add_argument("--section", default="Input")
    parser.add_argument("--sheet", default="single-hit")
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument("--title", default="Kang24 single-hit")
    parser.add_argument("--delta-activity", action="store_true")
    args = parser.parse_args()

    observed = read_observed(args.sheet)
    predicted = predictions_to_long(run_unfold_rate(args.xml, args.section))
    merged = observed.merge(predicted, on="xml_name", how="inner")
    if args.delta_activity:
        merged = add_delta_activity(merged)
        xcol, ycol = "observed_delta", "predicted_delta"
        xlabel, ylabel = "Observed Dactivity, log2(variant / WT)", "Predicted Dactivity, log2(variant / WT)"
        profile_y = "Mean Dactivity"
    else:
        xcol, ycol = "expression level", "predicted"
        xlabel, ylabel = "Observed expression", "Predicted expression"
        profile_y = "Mean expression"
    args.prefix.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.prefix.with_name(args.prefix.name + "_observed_predicted.csv"), index=False)
    make_scatter(merged, args.prefix.with_name(args.prefix.name + "_scatter.svg"), args.title, xcol, ycol, xlabel, ylabel)
    make_profile(merged, args.prefix.with_name(args.prefix.name + "_profile.svg"), args.title, xcol, ycol, profile_y)
    print(f"merged rows: {len(merged)}")
    print(f"pearson r: {merged[xcol].corr(merged[ycol]):.6f}")
    print(f"wrote {args.prefix.parent}")


if __name__ == "__main__":
    main()
