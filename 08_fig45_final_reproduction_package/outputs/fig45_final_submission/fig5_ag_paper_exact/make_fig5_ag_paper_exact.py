#!/usr/bin/env python3
"""Create paper-style Figure 5A/G panels from reproduced A/T substitution outputs."""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "outputs" / "fig45_competition_off_repeats" / "4tf_seed0" / "fig5" / "kang24_fig5_AT_substitution_table.csv"
OUT = ROOT / "outputs" / "fig45_final_submission" / "fig5_ag_paper_exact"

BLUE = "#001EFF"
RED = "#FF1010"
GRID = "#E7E7E7"
MINOR_GRID = "#F2F2F2"
AXIS = "#555555"
CRE_FILL = "#FFB6BD"
CRYPTIC_FILL = "#8F91FF"

REGIONS = [
    ("CRE1", 11, 18, CRE_FILL),
    ("CRE2", 37, 41, CRE_FILL),
    ("CRE3", 47, 51, CRE_FILL),
    ("cryptic", 63, 67, CRYPTIC_FILL),
    ("CRE4", 69, 76, CRE_FILL),
]
MUT_LABELS = {
    "A": [(11, "11A"), (14, "14A"), (17, "17A"), (69, "69A"), (72, "72A"), (75, "75A")],
    "T": [(12, "12T"), (15, "15T"), (18, "18T"), (70, "70T"), (73, "73T"), (76, "76T")],
}


def esc(s: object) -> str:
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def scale(v: float, a: float, b: float, c: float, d: float) -> float:
    if a == b:
        return (c + d) / 2
    return c + (v - a) / (b - a) * (d - c)


def text(x: float, y: float, s: str, size: int = 15, weight: str = "400", anchor: str = "start", fill: str = "#111111") -> str:
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" fill="{fill}">{esc(s)}</text>'
    )


def vtext(x: float, y: float, s: str, size: int = 17) -> str:
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-family="Arial, Helvetica, sans-serif" font-size="{size}" '
        f'text-anchor="middle" fill="#111111" transform="rotate(-90 {x:.2f} {y:.2f})">{esc(s)}</text>'
    )


def line(x1: float, y1: float, x2: float, y2: float, color: str = AXIS, width: float = 1.0) -> str:
    return f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{color}" stroke-width="{width}"/>'


def rect(x: float, y: float, w: float, h: float, fill: str, stroke: str = "none", width: float = 1.0, opacity: float = 1.0) -> str:
    return f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" fill="{fill}" stroke="{stroke}" stroke-width="{width}" opacity="{opacity}"/>'


def write_svg(path: Path, width: int, height: int, body: list[str], scale_factor: int = 1) -> None:
    path.write_text(
        "\n".join(
            [
                f'<svg xmlns="http://www.w3.org/2000/svg" width="{width * scale_factor}" height="{height * scale_factor}" viewBox="0 0 {width} {height}">',
                '<rect width="100%" height="100%" fill="#ffffff"/>',
                *body,
                "</svg>",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def stats(df: pd.DataFrame, base: str) -> tuple[float, float, float]:
    sub = df[df["mut_base"].eq(base)].dropna(subset=["observed_delta", "predicted_delta"])
    r = float(sub["observed_delta"].corr(sub["predicted_delta"]))
    r2 = r * r
    rms = math.sqrt(float(((sub["observed_delta"] - sub["predicted_delta"]) ** 2).mean()))
    return r, r2, rms


def draw_axis(body: list[str], x: float, y: float, w: float, h: float, show_x: bool) -> None:
    xlo, xhi = 0, 86
    ylo, yhi = -2.2, 2.2
    body.append(rect(x, y, w, h, "#FFFFFF", "#777777", 0.9))
    for label, lo, hi, color in REGIONS:
        rx = scale(lo, xlo, xhi, x, x + w)
        rw = scale(hi, xlo, xhi, x, x + w) - rx
        body.append(rect(rx, y, rw, h, color, "none", 1, 0.64))
    for tick in range(0, 90, 10):
        tx = scale(tick, xlo, xhi, x, x + w)
        body.append(line(tx, y, tx, y + h, GRID if tick % 20 == 0 else MINOR_GRID, 0.8))
        if show_x:
            body.append(text(tx, y + h + 21, f"{tick}", 13, anchor="middle", fill="#666666"))
    for tick in [-2, -1, 0, 1, 2]:
        ty = scale(tick, ylo, yhi, y + h, y)
        body.append(line(x, ty, x + w, ty, GRID, 0.9))
        body.append(text(x - 7, ty + 4, f"{tick:g}", 13, anchor="end", fill="#555555"))
    zero = scale(0, ylo, yhi, y + h, y)
    body.append(line(x, zero, x + w, zero, "#222222", 0.9))


def draw_bars(body: list[str], sub: pd.DataFrame, value_col: str, x: float, y: float, w: float, h: float, color: str) -> None:
    xlo, xhi = 0, 86
    ylo, yhi = -2.2, 2.2
    zero = scale(0, ylo, yhi, y + h, y)
    bar_w = max(2.0, w / 86 * 0.55)
    for row in sub.itertuples(index=False):
        pos = float(row.position)
        value = float(getattr(row, value_col))
        px = scale(pos, xlo, xhi, x, x + w)
        py = scale(value, ylo, yhi, y + h, y)
        body.append(rect(px - bar_w / 2, min(zero, py), bar_w, abs(zero - py), color, "none"))


def draw_region_labels(body: list[str], x: float, y: float, w: float) -> None:
    for label, lo, hi, _ in REGIONS:
        lx = scale((lo + hi) / 2, 0, 86, x, x + w)
        body.append(text(lx, y - 5, label, 14, anchor="middle"))


def draw_mutation_brackets(body: list[str], base: str, x: float, y: float, w: float, h: float) -> None:
    y_base = y + h - 42
    offsets = {
        "A": {11: -10, 14: 0, 17: 10, 69: -10, 72: 0, 75: 10},
        "T": {12: -10, 15: 0, 18: 10, 70: -10, 73: 0, 76: 10},
    }
    for pos, label in MUT_LABELS[base]:
        px = scale(pos, 0, 86, x, x + w)
        body.append(line(px, y_base - 10, px, y_base + 7, "#333333", 0.8))
        label_x = px + offsets[base][pos]
        body.append(line(px, y_base + 7, label_x, y_base + 20, "#333333", 0.8))
        body.append(text(label_x, y_base + 33, label, 11, anchor="middle"))


def panel(body: list[str], df: pd.DataFrame, base: str, letter: str, x: float, y: float, w: float, panel_h: float) -> None:
    gap = 14
    top_h = panel_h
    bottom_h = panel_h
    sub = df[df["mut_base"].eq(base)].dropna(subset=["position"]).copy()
    body.append(text(x - 67, y - 2, letter, 18, "700"))
    draw_axis(body, x, y, w, top_h, False)
    draw_axis(body, x, y + top_h + gap, w, bottom_h, True)
    draw_region_labels(body, x, y, w)
    draw_bars(body, sub, "observed_delta", x, y, w, top_h, BLUE)
    draw_bars(body, sub, "predicted_delta", x, y + top_h + gap, w, bottom_h, RED)
    body.append(text(x + 24, y + 24, f"→{base}", 22))
    body.append(text(x + 24, y + top_h + gap + 24, f"→{base}", 22))
    r, r2, rms = stats(df, base)
    body.append(text(x + w - 45, y + top_h + gap + 25, f"R={r:.2f} R²={r2:.2f} rms={rms:.2g}", 15, anchor="end"))
    draw_mutation_brackets(body, base, x, y + top_h + gap, w, bottom_h)
    body.append(vtext(x - 34, y + top_h / 2, "Δactivity(log₂)", 16))
    body.append(vtext(x - 34, y + top_h + gap + bottom_h / 2, "Δactivity(log₂)", 16))
    body.append(text(x + w / 2, y + top_h + gap + bottom_h + 50, "position", 15, anchor="middle"))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(DATA)
    width, height = 1280, 390
    body: list[str] = []
    panel(body, df, "A", "A", 95, 40, 500, 122)
    panel(body, df, "T", "G", 735, 40, 500, 122)
    svg = OUT / "fig5_ag_paper_exact.svg"
    svg_4x = OUT / "fig5_ag_paper_exact_4x.svg"
    write_svg(svg, width, height, body)
    write_svg(svg_4x, width, height, body, scale_factor=4)
    rows = []
    for base, letter in [("A", "A"), ("T", "G")]:
        r, r2, rms = stats(df, base)
        rows.append({"panel": letter, "base": base, "n": int(df["mut_base"].eq(base).sum()), "pearson_r": r, "r_squared": r2, "rms": rms})
    pd.DataFrame(rows).to_csv(OUT / "fig5_ag_paper_exact_metrics.csv", index=False)
    print(svg)


if __name__ == "__main__":
    main()
