#!/usr/bin/env python3
"""Create a paper-style Figure 4J panel from reproduced multi-hit outputs."""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "fig45_final_submission" / "fig4j_paper_exact"
DATA = ROOT / "outputs" / "fig45_competition_off_repeats" / "4tf_seed0" / "fig4j" / "kang24_fig4j_4tf_multihit_observed_predicted.csv"

BLUE = "#1f4fb2"
PANEL_BG = "#EBEBEB"
GRID = "#FFFFFF"
AXIS = "#555555"


def esc(s: object) -> str:
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def scale(v: float, src_min: float, src_max: float, dst_min: float, dst_max: float) -> float:
    if src_max == src_min:
        return (dst_min + dst_max) / 2
    return dst_min + (v - src_min) / (src_max - src_min) * (dst_max - dst_min)


def text(x: float, y: float, s: str, size: int = 18, weight: str = "400", anchor: str = "start", fill: str = "#111111") -> str:
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" fill="{fill}">{esc(s)}</text>'
    )


def vtext(x: float, y: float, s: str, size: int = 25, weight: str = "400") -> str:
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" text-anchor="middle" fill="#111111" '
        f'transform="rotate(-90 {x:.2f} {y:.2f})">{esc(s)}</text>'
    )


def line(x1: float, y1: float, x2: float, y2: float, color: str = AXIS, width: float = 1.0) -> str:
    return f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{color}" stroke-width="{width}"/>'


def rect(x: float, y: float, w: float, h: float, fill: str, stroke: str = "none", width: float = 1.0) -> str:
    return f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" fill="{fill}" stroke="{stroke}" stroke-width="{width}"/>'


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


def make_panel(df: pd.DataFrame) -> list[str]:
    # Match the published Fig. 4J axis window: x around -4..1, y around -2.7..0.9.
    xlo, xhi = -4.1, 1.1
    ylo, yhi = -2.8, 0.9
    x, y, w, h = 100, 72, 285, 285
    body: list[str] = []
    body.append(text(34, 48, "J", 20, "700"))
    body.append(rect(x, y, w, h, PANEL_BG))

    x_major = [-4, -2, 0]
    x_minor = [-3, -1, 1]
    y_major = [-2, -1, 0]
    y_minor = [-2.5, -1.5, -0.5, 0.5]
    for tick in x_minor:
        tx = scale(tick, xlo, xhi, x, x + w)
        body.append(line(tx, y, tx, y + h, GRID, 0.85))
    for tick in y_minor:
        ty = scale(tick, ylo, yhi, y + h, y)
        body.append(line(x, ty, x + w, ty, GRID, 0.85))
    for tick in x_major:
        tx = scale(tick, xlo, xhi, x, x + w)
        body.append(line(tx, y, tx, y + h, GRID, 1.45))
        body.append(line(tx, y + h, tx, y + h + 6, AXIS, 1.1))
        body.append(text(tx, y + h + 28, f"{tick:g}", 18, anchor="middle", fill="#4c4c4c"))
    for tick in y_major:
        ty = scale(tick, ylo, yhi, y + h, y)
        body.append(line(x, ty, x + w, ty, GRID, 1.45))
        body.append(line(x - 6, ty, x, ty, AXIS, 1.1))
        body.append(text(x - 10, ty + 6, f"{tick:g}", 18, anchor="end", fill="#4c4c4c"))

    body.append(line(x, y + h, x + w, y + h, AXIS, 1.1))
    body.append(line(x, y, x, y + h, AXIS, 1.1))
    for row in df.itertuples(index=False):
        px = scale(float(row.observed_delta), xlo, xhi, x, x + w)
        py = scale(float(row.predicted_delta), ylo, yhi, y + h, y)
        if x - 5 <= px <= x + w + 5 and y - 5 <= py <= y + h + 5:
            body.append(f'<circle cx="{px:.2f}" cy="{py:.2f}" r="1.45" fill="{BLUE}" opacity="0.86"/>')

    r_value = float(df.observed_delta.corr(df.predicted_delta))
    body.append(text(x + 8, y + 28, f"R={r_value:.2f}", 19))
    body.append(vtext(45, y + h / 2, "Multi-hit prediction", 25))
    body.append(text(x + w / 2, y + h + 65, "MPRA multi-hit data", 25, anchor="middle"))
    return body


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(DATA)
    body = make_panel(df)
    width, height = 420, 430
    svg = OUT / "fig4j_our_model_paper_exact.svg"
    svg_4x = OUT / "fig4j_our_model_paper_exact_4x.svg"
    write_svg(svg, width, height, body)
    write_svg(svg_4x, width, height, body, scale_factor=4)
    summary = pd.DataFrame(
        [
            {
                "panel": "J",
                "n": len(df),
                "pearson_r": df.observed_delta.corr(df.predicted_delta),
                "rmse": math.sqrt(float(((df.observed_delta - df.predicted_delta) ** 2).mean())),
                "x_min": df.observed_delta.min(),
                "x_max": df.observed_delta.max(),
                "y_min": df.predicted_delta.min(),
                "y_max": df.predicted_delta.max(),
            }
        ]
    )
    summary.to_csv(OUT / "fig4j_our_model_paper_exact_metrics.csv", index=False)
    print(svg)


if __name__ == "__main__":
    main()
