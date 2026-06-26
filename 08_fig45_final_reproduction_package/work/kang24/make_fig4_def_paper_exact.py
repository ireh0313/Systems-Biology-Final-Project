#!/usr/bin/env python3
"""Create a paper-style Figure 4D-F panel from reproduced model outputs."""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
FINAL = ROOT / "outputs" / "fig45_final_submission"
OUT = ROOT / "outputs" / "fig45_final_submission" / "fig4_def_paper_exact"

CASE_FILES = [
    ("D", "Case 1", FINAL / "fig4_reverse_rearrange" / "kang24_fig4_reverse_observed_predicted.csv"),
    ("E", "Case 2", FINAL / "fig4_reverse_rearrange" / "kang24_fig4_rearrange_observed_predicted.csv"),
    ("F", "Case 3", FINAL / "fig4_reverse_rearrange" / "kang24_fig4_reverse_and_rearrange_observed_predicted.csv"),
]
BASELINE = ROOT / "outputs" / "fig45_competition_off_repeats" / "4tf_seed0" / "kang24_4tf_seed0_singlehit_delta_observed_predicted.csv"

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


def vtext(x: float, y: float, s: str, size: int = 21, weight: str = "400") -> str:
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
    canvas_w = width * scale_factor
    canvas_h = height * scale_factor
    path.write_text(
        "\n".join(
            [
                f'<svg xmlns="http://www.w3.org/2000/svg" width="{canvas_w}" height="{canvas_h}" viewBox="0 0 {width} {height}">',
                '<rect width="100%" height="100%" fill="#ffffff"/>',
                *body,
                "</svg>",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def points(df: pd.DataFrame, x: float, y: float, w: float, h: float, r: float, alpha: float = 0.9) -> list[str]:
    xlo, xhi = -2.1, 1.1
    ylo, yhi = -2.1, 1.1
    out = []
    for row in df.itertuples(index=False):
        px = scale(float(row.observed_delta), xlo, xhi, x, x + w)
        py = scale(float(row.predicted_delta), ylo, yhi, y + h, y)
        out.append(f'<circle cx="{px:.2f}" cy="{py:.2f}" r="{r:.2f}" fill="{BLUE}" opacity="{alpha:.2f}"/>')
    return out


def scatter_panel(letter: str, df: pd.DataFrame, baseline: pd.DataFrame, x: float, y: float, w: float, h: float, show_xlabel: bool) -> list[str]:
    xlo, xhi = -2.1, 1.1
    ylo, yhi = -2.1, 1.1
    body: list[str] = []
    body.append(text(x - 64, y - 16, letter, 18, "700"))
    body.append(rect(x, y, w, h, PANEL_BG, "none"))

    major_ticks = [-2, -1, 0, 1]
    minor_ticks = [-1.5, -0.5, 0.5]
    for tick in minor_ticks:
        tx = scale(tick, xlo, xhi, x, x + w)
        ty = scale(tick, ylo, yhi, y + h, y)
        body.append(line(tx, y, tx, y + h, GRID, 0.75))
        body.append(line(x, ty, x + w, ty, GRID, 0.75))
    for tick in major_ticks:
        tx = scale(tick, xlo, xhi, x, x + w)
        ty = scale(tick, ylo, yhi, y + h, y)
        body.append(line(tx, y, tx, y + h, GRID, 1.4))
        body.append(line(x, ty, x + w, ty, GRID, 1.4))
        body.append(line(tx, y + h, tx, y + h + 5, AXIS, 1.0))
        body.append(line(x - 5, ty, x, ty, AXIS, 1.0))
        body.append(text(tx, y + h + 25, f"{tick:.0f}", 16, anchor="middle", fill="#4c4c4c"))
        body.append(text(x - 10, ty + 6, f"{tick:.1f}", 16, anchor="end", fill="#4c4c4c"))
    body.append(line(x, y + h, x + w, y + h, AXIS, 1.0))
    body.append(line(x, y, x, y + h, AXIS, 1.0))
    body += points(df, x, y, w, h, 1.25, 0.82)

    inset_x = x + 13
    inset_y = y + 13
    inset_w = 87
    inset_h = 87
    body.append(rect(inset_x, inset_y, inset_w, inset_h, "#F2F2F2", "#A0A0A0", 0.8))
    for tick in major_ticks:
        tx = scale(tick, xlo, xhi, inset_x, inset_x + inset_w)
        ty = scale(tick, ylo, yhi, inset_y + inset_h, inset_y)
        body.append(line(tx, inset_y, tx, inset_y + inset_h, GRID, 0.55))
        body.append(line(inset_x, ty, inset_x + inset_w, ty, GRID, 0.55))
    body.append(line(inset_x, inset_y + inset_h, inset_x + inset_w, inset_y, BLUE, 0.9))
    body += points(baseline, inset_x, inset_y, inset_w, inset_h, 0.65, 0.78)
    if letter == "D":
        body.append(text(inset_x + 8, inset_y + inset_h - 10, f"R={baseline.observed_delta.corr(baseline.predicted_delta):.2f}", 15))

    r_value = float(df.observed_delta.corr(df.predicted_delta))
    body.append(text(x + w - 8, y + h - 19, f"R={r_value:.3f}", 17, anchor="end"))
    body.append(vtext(x - 53, y + h / 2, "prediction", 21))
    if show_xlabel:
        body.append(text(x + w / 2, y + h + 57, "MPRA single-hit data", 21, anchor="middle"))
    return body


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    baseline = pd.read_csv(BASELINE)
    cases = [(letter, pd.read_csv(path)) for letter, _, path in CASE_FILES]

    width, height = 405, 930
    plot_x, plot_w, plot_h = 118, 235, 205
    ys = [82, 352, 622]
    body = [text(plot_x + plot_w / 2, 42, "Our model", 24, "400", anchor="middle")]
    for (letter, df), y in zip(cases, ys):
        body += scatter_panel(letter, df, baseline, plot_x, y, plot_w, plot_h, show_xlabel=(letter == "F"))

    out_svg = OUT / "fig4_def_our_model_paper_exact.svg"
    out_svg_4x = OUT / "fig4_def_our_model_paper_exact_4x.svg"
    write_svg(out_svg, width, height, body)
    write_svg(out_svg_4x, width, height, body, scale_factor=4)

    summary = []
    for letter, df in cases:
        summary.append(
            {
                "panel": letter,
                "n": len(df),
                "pearson_r": df.observed_delta.corr(df.predicted_delta),
                "rmse": math.sqrt(float(((df.observed_delta - df.predicted_delta) ** 2).mean())),
            }
        )
    pd.DataFrame(summary).to_csv(OUT / "fig4_def_our_model_paper_exact_metrics.csv", index=False)
    print(out_svg)


if __name__ == "__main__":
    main()
