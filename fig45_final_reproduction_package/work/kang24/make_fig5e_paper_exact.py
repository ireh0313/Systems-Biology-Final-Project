#!/usr/bin/env python3
"""Create a paper-style focused reproduction of Kang et al. Figure 5E."""

from __future__ import annotations

from pathlib import Path
import math
import sys

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
KANG = ROOT / "work" / "kang24"
sys.path.insert(0, str(KANG))
import make_fig5_analysis as fig5  # noqa: E402


XML = KANG / "fig45_competition_off_repeats" / "kang24_4tf_seed0.xml"
OUT = ROOT / "outputs" / "fig45_final_submission" / "fig5e_paper_exact"

TF_ORDER = ["CREB1", "CREM", "ATF1", "ATF7"]
TF_COLORS = {
    "CREB1": "#FF9900",
    "CREM": "#49A000",
    "ATF1": "#20208F",
    "ATF7": "#B442FF",
}


def esc(s: object) -> str:
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def scale(v: float, a: float, b: float, c: float, d: float) -> float:
    if a == b:
        return (c + d) / 2
    return c + (v - a) / (b - a) * (d - c)


def text(
    x: float,
    y: float,
    s: str,
    size: int = 16,
    weight: str = "400",
    anchor: str = "start",
    fill: str = "#111111",
    extra: str = "",
) -> str:
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" fill="{fill}" {extra}>{esc(s)}</text>'
    )


def line(x1: float, y1: float, x2: float, y2: float, color: str = "#222222", width: float = 1.0) -> str:
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


def rotated_label(cx: float, cy: float, label: str, size: int = 14) -> str:
    return text(cx, cy, label, size, anchor="end", extra=f'transform="rotate(-45 {cx:.2f} {cy:.2f})"')


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    coefs = fig5.coefficients(XML)
    rows = []
    for tf in TF_ORDER:
        coef = float(coefs[tf])
        rows.append({"tf": tf, "coefficient": coef, "log10_coefficient": math.log10(coef)})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "fig5e_activation_coefficients.csv", index=False)

    width, height = 300, 340
    left, right, top, bottom = 80, 24, 20, 108
    plot_w = width - left - right
    plot_h = height - top - bottom
    ymin, ymax = -2.2, 1.35
    zero_y = scale(0, ymin, ymax, top + plot_h, top)
    body: list[str] = []
    body.append(text(15, 22, "E", 20, "700"))

    # Axes and light grid match the compact panel style of the paper.
    axis_x = left
    axis_bottom = top + plot_h
    body.append(line(axis_x, top, axis_x, axis_bottom, "#222222", 1.1))
    body.append(line(axis_x, axis_bottom, left + plot_w + 4, axis_bottom, "#222222", 1.1))
    for tick in [0.8, -0.2, -1.2, -2.2]:
        ty = scale(tick, ymin, ymax, axis_bottom, top)
        body.append(line(axis_x - 5, ty, axis_x, ty, "#222222", 1.0))
        body.append(text(axis_x - 9, ty + 4, f"{tick:.1f}", 15, anchor="end"))

    # The paper shows bars from the bottom baseline rather than from zero.
    bar_w = 29
    gap = 16
    start_x = left + 16
    baseline = axis_bottom
    for i, row in enumerate(df.itertuples(index=False)):
        tf = str(row.tf)
        value = float(row.log10_coefficient)
        bx = start_x + i * (bar_w + gap)
        clipped = min(ymax, max(ymin, value))
        by = scale(clipped, ymin, ymax, axis_bottom, top)
        if value < ymin:
            by = axis_bottom - 10
        body.append(rect(bx, by, bar_w, baseline - by, TF_COLORS[tf], "none"))
        tick_x = bx + bar_w / 2
        body.append(line(tick_x, axis_bottom, tick_x, axis_bottom + 5, "#222222", 1.0))
        body.append(rotated_label(bx + bar_w * 0.62, height - 44, tf, 15))

    body.append(
        text(
            18,
            top + plot_h / 2,
            "log10(coef)",
            16,
            anchor="middle",
            extra=f'transform="rotate(-90 18.00 {top + plot_h / 2:.2f})"',
        )
    )
    write_svg(OUT / "fig5e_paper_exact.svg", width, height, body)
    write_svg(OUT / "fig5e_paper_exact_4x.svg", width, height, body, scale_factor=4)
    print(OUT / "fig5e_paper_exact.svg")


if __name__ == "__main__":
    main()
