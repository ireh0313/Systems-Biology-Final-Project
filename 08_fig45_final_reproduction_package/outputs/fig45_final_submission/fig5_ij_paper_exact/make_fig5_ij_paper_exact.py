#!/usr/bin/env python3
"""Create paper-style Figure 5I/J T-substitution zoom panels."""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
KANG = ROOT / "work" / "kang24"
sys.path.insert(0, str(KANG))
import make_fig5_analysis as fig5  # noqa: E402


XML = KANG / "fig45_competition_off_repeats" / "kang24_4tf_seed0.xml"
AT_TABLE = ROOT / "outputs" / "fig45_competition_off_repeats" / "4tf_seed0" / "fig5" / "kang24_fig5_AT_substitution_table.csv"
OUT = ROOT / "outputs" / "fig45_final_submission" / "fig5_ij_paper_exact"

TF_ORDER = ["ATF7", "ATF1", "CREM", "CREB1"]
TF_COLORS = {
    "ATF7": "#B442FF",
    "ATF1": "#20208F",
    "CREM": "#49A000",
    "CREB1": "#FF9900",
}
REGION_COLOR = "#FF1F2D"
ALT_REGION_COLOR = "#1920D7"
MUT_COLOR = "#FFF200"

PANELS = [
    ("I", "CRE1", "12T", "scanmut_single_pos_12_T", 12, 8, 22, 20, 45),
    ("I", "CRE1", "15T", "scanmut_single_pos_15_T", 15, 8, 22, 20, 45),
    ("I", "CRE1", "18T", "scanmut_single_pos_18_T", 18, 8, 22, 20, 45),
    ("J", "CRE4", "70T", "scanmut_single_pos_70_T", 70, 66, 80, 78, 102),
    ("J", "CRE4", "73T", "scanmut_single_pos_73_T", 73, 66, 80, 78, 102),
    ("J", "CRE4", "76T", "scanmut_single_pos_76_T", 76, 66, 80, 78, 102),
]


def esc(s: object) -> str:
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def scale(v: float, a: float, b: float, c: float, d: float) -> float:
    if a == b:
        return (c + d) / 2
    return c + (v - a) / (b - a) * (d - c)


def text(x: float, y: float, s: str, size: int = 16, weight: str = "400", anchor: str = "start", fill: str = "#111111") -> str:
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
    return 0.14 + 0.66 * min(1.0, max(0.0, value / max_occ if max_occ else 0))


def sequence_fragment(at: pd.DataFrame, gene: str, scan_lo: int, scan_hi: int) -> str:
    row = at[at["xml_name"].eq(gene)].iloc[0]
    seq = str(row["sequence"])
    return seq[20 + scan_lo : 20 + scan_hi]


def draw_sequence(body: list[str], frag: str, mut_pos: int, scan_lo: int, x: float, y: float, w: float) -> None:
    char_w = w / len(frag)
    mut_idx = int(mut_pos - scan_lo)
    red_start = max(0, (11 if scan_lo < 40 else 69) - scan_lo)
    red_end = min(len(frag), (19 if scan_lo < 40 else 77) - scan_lo)
    body.append(rect(x + red_start * char_w - 1, y - 15, (red_end - red_start) * char_w + 2, 20, "none", REGION_COLOR, 1.5))
    body.append(rect(x + mut_idx * char_w, y - 15, char_w, 20, MUT_COLOR, "none", opacity=0.85))
    for i, ch in enumerate(frag):
        body.append(text(x + (i + 0.5) * char_w, y, ch, 16, anchor="middle"))


def draw_panel(
    body: list[str],
    sites: pd.DataFrame,
    at: pd.DataFrame,
    label: str,
    gene: str,
    mut_pos: int,
    scan_lo: int,
    scan_hi: int,
    xlo: int,
    xhi: int,
    x: float,
    y: float,
    w: float,
    h: float,
    show_letter: str | None = None,
) -> None:
    body.append(rect(x, y, w, h, "#FFFFFF", "#222222", 1.2))
    if show_letter:
        body.append(text(x - 88, y - 10, show_letter, 20, "700"))
    body.append(text(x + w - 8, y + 24, label, 22, anchor="end"))
    row_gap = 28
    box_h = 24
    max_occ = float(sites["occupancy"].max())
    for i, tf in enumerate(TF_ORDER):
        yy = y + 27 + i * row_gap
        body.append(text(x - 12, yy + 6, tf, 17, anchor="end"))
        sub = sites[(sites["gene"].eq(gene)) & (sites["tf"].eq(tf))]
        for row in sub.itertuples(index=False):
            if float(row.seq_end) < xlo or float(row.seq_start) > xhi:
                continue
            sx = scale(max(float(row.seq_start), xlo), xlo, xhi, x, x + w)
            ex = scale(min(float(row.seq_end), xhi), xlo, xhi, x, x + w)
            bw = max(3.0, ex - sx)
            opacity = occupancy_alpha(float(row.occupancy), max_occ)
            color = TF_COLORS[tf]
            body.append(rect(sx, yy - box_h / 2, bw, box_h, color, "none", opacity=opacity))
            body.append(rect(sx, yy - box_h / 2, bw, box_h, "none", color, 2.0))

    base_y = y + h - 16
    body.append(line(x, base_y, x + w, base_y, "#333333", 1.0))
    cre_lo, cre_hi = ((11, 19) if xlo < 50 else (69, 77))
    if xlo >= 50:
        bx = scale(62 + 20, xlo, xhi, x, x + w)
        bw = scale(65 + 20, xlo, xhi, x, x + w) - bx
        body.append(rect(bx, base_y - 8, bw, 16, ALT_REGION_COLOR, "none", opacity=0.95))
    rx = scale(cre_lo + 20, xlo, xhi, x, x + w)
    rw = scale(cre_hi + 20, xlo, xhi, x, x + w) - rx
    body.append(rect(rx, base_y - 8, rw, 16, REGION_COLOR, "none", opacity=0.95))
    mx = scale(mut_pos + 20, xlo, xhi, x, x + w)
    body.append(rect(mx - 2.2, base_y - 10, 4.4, 20, MUT_COLOR, "none", opacity=0.95))
    seq_y = y + h + 24
    body.append(line(rx, base_y + 8, x + w * 0.30, seq_y - 17, "#222222", 1.0))
    body.append(line(rx + rw, base_y + 8, x + w * 0.70, seq_y - 17, "#222222", 1.0))
    frag = sequence_fragment(at, gene, scan_lo, scan_hi)
    draw_sequence(body, frag, mut_pos, scan_lo, x - 4, seq_y, w + 8)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    at = pd.read_csv(AT_TABLE)
    genes = [p[3] for p in PANELS]
    frames = []
    for gene in genes:
        s = fig5.sites_for_gene(XML, gene)
        s["gene"] = gene
        frames.append(s)
    sites = pd.concat(frames, ignore_index=True)
    sites.to_csv(OUT / "fig5_ij_variant_tfbs_occupancy_sites.csv", index=False)

    width, height = 620, 680
    panel_w, panel_h = 170, 165
    x_left, x_right = 115, 405
    ys = [32, 244, 456]
    body: list[str] = []
    for idx, spec in enumerate(PANELS[:3]):
        letter, _, label, gene, mut_pos, scan_lo, scan_hi, xlo, xhi = spec
        draw_panel(body, sites, at, label, gene, mut_pos, scan_lo, scan_hi, xlo, xhi, x_left, ys[idx], panel_w, panel_h, show_letter=(letter if idx == 0 else None))
    for idx, spec in enumerate(PANELS[3:]):
        letter, _, label, gene, mut_pos, scan_lo, scan_hi, xlo, xhi = spec
        draw_panel(body, sites, at, label, gene, mut_pos, scan_lo, scan_hi, xlo, xhi, x_right, ys[idx], panel_w, panel_h, show_letter=(letter if idx == 0 else None))

    svg = OUT / "fig5_ij_paper_exact.svg"
    svg_4x = OUT / "fig5_ij_paper_exact_4x.svg"
    write_svg(svg, width, height, body)
    write_svg(svg_4x, width, height, body, scale_factor=4)
    summary = sites.groupby(["gene", "tf"], as_index=False).agg(n_sites=("occupancy", "size"), max_occupancy=("occupancy", "max"))
    summary.to_csv(OUT / "fig5_ij_paper_exact_metrics.csv", index=False)
    print(svg)


if __name__ == "__main__":
    main()
