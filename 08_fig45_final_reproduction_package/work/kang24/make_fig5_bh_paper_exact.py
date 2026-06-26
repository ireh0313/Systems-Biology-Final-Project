#!/usr/bin/env python3
"""Create paper-style Figure 5B/H TFBS occupancy maps from reproduced WT output."""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
KANG = ROOT / "work" / "kang24"
sys.path.insert(0, str(KANG))
import make_fig5_analysis as fig5  # noqa: E402


XML = KANG / "fig45_competition_off_repeats" / "kang24_4tf_seed0.xml"
OUT = ROOT / "outputs" / "fig45_final_submission" / "fig5_bh_paper_exact"

TF_ORDER = ["ATF7", "ATF1", "CREM", "CREB1"]
TF_COLORS = {
    "ATF7": "#B442FF",
    "ATF1": "#20208F",
    "CREM": "#49A000",
    "CREB1": "#FF9900",
}
REGIONS = {
    "CRE1": (11, 18, "#FF1F2D"),
    "CRE2": (37, 41, "#FF1F2D"),
    "CRE3": (47, 51, "#FF1F2D"),
    "cryptic": (63, 67, "#3030DD"),
    "CRE4": (69, 76, "#FF1F2D"),
}


def esc(s: object) -> str:
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def scale(v: float, a: float, b: float, c: float, d: float) -> float:
    if a == b:
        return (c + d) / 2
    return c + (v - a) / (b - a) * (d - c)


def text(x: float, y: float, s: str, size: int = 17, weight: str = "400", anchor: str = "start", fill: str = "#111111") -> str:
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" fill="{fill}">{esc(s)}</text>'
    )


def line(x1: float, y1: float, x2: float, y2: float, color: str = "#222222", width: float = 1.0) -> str:
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


def occupancy_alpha(value: float, max_occ: float) -> float:
    if max_occ <= 0:
        return 0.18
    return 0.16 + 0.64 * min(1.0, max(0.0, value / max_occ))


def draw_panel(body: list[str], sites: pd.DataFrame, letter: str, x: float, y: float, w: float, h: float) -> None:
    xlo, xhi = 0, 110
    inner_x = x + 34
    inner_w = w - 48
    row_gap = 29
    box_h = 24
    body.append(text(x - 54, y - 1, letter, 19, "700"))
    body.append(rect(inner_x, y, inner_w, h, "#FFFFFF", "#222222", 1.2))

    max_occ = float(sites["occupancy"].max())
    row_centers: dict[str, float] = {}
    for i, tf in enumerate(TF_ORDER):
        yy = y + 24 + i * row_gap
        row_centers[tf] = yy
        body.append(text(x + 16, yy + 6, tf, 17, anchor="end"))
        tf_sites = sites[sites["tf"].eq(tf)].copy()
        for row in tf_sites.itertuples(index=False):
            if float(row.seq_end) < xlo or float(row.seq_start) > xhi:
                continue
            sx = scale(max(float(row.seq_start), xlo), xlo, xhi, inner_x, inner_x + inner_w)
            ex = scale(min(float(row.seq_end), xhi), xlo, xhi, inner_x, inner_x + inner_w)
            opacity = occupancy_alpha(float(row.occupancy), max_occ)
            color = TF_COLORS[tf]
            bw = max(2.5, ex - sx)
            body.append(rect(sx, yy - box_h / 2, bw, box_h, color, "none", 1.0, opacity))
            body.append(rect(sx, yy - box_h / 2, bw, box_h, "none", color, 2.0, 1.0))

    base_y = y + h - 18
    body.append(line(inner_x, base_y, inner_x + inner_w, base_y, "#333333", 1.2))
    body.append(rect(inner_x, base_y - 5, 22, 10, "#C9C9C9", "none", opacity=0.9))
    body.append(rect(inner_x + inner_w - 22, base_y - 5, 22, 10, "#C9C9C9", "none", opacity=0.9))
    for name, (lo, hi, color) in REGIONS.items():
        sx = scale(lo + 20, xlo, xhi, inner_x, inner_x + inner_w)
        ex = scale(hi + 20, xlo, xhi, inner_x, inner_x + inner_w)
        body.append(rect(sx, base_y - 8, max(2, ex - sx), 16, color, "none", opacity=0.95))


def draw_fraction_legend(body: list[str], x: float, y: float) -> None:
    body.append(text(x, y, "High F", 16))
    grad_id = "frac_grad"
    body.append(
        f'<defs><linearGradient id="{grad_id}" x1="0" x2="0" y1="0" y2="1">'
        '<stop offset="0%" stop-color="#111111" stop-opacity="1"/>'
        '<stop offset="100%" stop-color="#111111" stop-opacity="0.08"/>'
        '</linearGradient></defs>'
    )
    body.append(f'<polygon points="{x+6:.2f},{y+20:.2f} {x+46:.2f},{y+20:.2f} {x+6:.2f},{y+120:.2f}" fill="url(#{grad_id})" stroke="#222222" stroke-width="0.8"/>')
    body.append(text(x, y + 142, "Low F", 16))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    sites = fig5.sites_for_gene(XML, "synCRE_Promega_0")
    sites.to_csv(OUT / "fig5_bh_wt_tfbs_occupancy_sites.csv", index=False)

    width, height = 1280, 270
    body: list[str] = []
    draw_panel(body, sites, "B", 90, 46, 500, 178)
    draw_panel(body, sites, "H", 690, 46, 500, 178)
    draw_fraction_legend(body, 1210, 45)

    svg = OUT / "fig5_bh_paper_exact.svg"
    svg_4x = OUT / "fig5_bh_paper_exact_4x.svg"
    write_svg(svg, width, height, body)
    write_svg(svg_4x, width, height, body, scale_factor=4)
    summary = sites.groupby("tf", as_index=False).agg(n_sites=("occupancy", "size"), mean_occupancy=("occupancy", "mean"), max_occupancy=("occupancy", "max"))
    summary.to_csv(OUT / "fig5_bh_paper_exact_metrics.csv", index=False)
    print(svg)


if __name__ == "__main__":
    main()
