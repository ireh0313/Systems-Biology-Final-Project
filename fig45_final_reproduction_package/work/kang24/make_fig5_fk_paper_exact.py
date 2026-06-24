#!/usr/bin/env python3
"""Create paper-style focused reproduction of Kang et al. Figure 5F/K."""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
KANG = ROOT / "work" / "kang24"
sys.path.insert(0, str(KANG))
import make_fig5_analysis as fig5  # noqa: E402


XML = KANG / "fig45_competition_off_repeats" / "kang24_4tf_seed0.xml"
OUT = ROOT / "outputs" / "fig45_final_submission" / "fig5_fk_paper_exact"

TF_STACK = ["CREB1", "CREM", "ATF1", "ATF7"]
TF_COLORS = {
    "CREB1": "#FF9900",
    "CREM": "#49A000",
    "ATF1": "#20208F",
    "ATF7": "#B442FF",
}

PANELS = {
    "F": {
        "base": "A",
        "left": [("WT", "synCRE_Promega_0"), ("11A", "scanmut_single_pos_11_A"), ("14A", "scanmut_single_pos_14_A"), ("17A", "scanmut_single_pos_17_A")],
        "right": [("WT", "synCRE_Promega_0"), ("69A", "scanmut_single_pos_69_A"), ("72A", "scanmut_single_pos_72_A"), ("75A", "scanmut_single_pos_75_A")],
    },
    "K": {
        "base": "T",
        "left": [("WT", "synCRE_Promega_0"), ("12T", "scanmut_single_pos_12_T"), ("15T", "scanmut_single_pos_15_T"), ("18T", "scanmut_single_pos_18_T")],
        "right": [("WT", "synCRE_Promega_0"), ("70T", "scanmut_single_pos_70_T"), ("73T", "scanmut_single_pos_73_T"), ("76T", "scanmut_single_pos_76_T")],
    },
}


def esc(s: object) -> str:
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def scale(v: float, a: float, b: float, c: float, d: float) -> float:
    if a == b:
        return (c + d) / 2
    return c + (v - a) / (b - a) * (d - c)


def text(x: float, y: float, s: str, size: int = 16, weight: str = "400", anchor: str = "start", fill: str = "#111111", extra: str = "") -> str:
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" fill="{fill}" {extra}>{esc(s)}</text>'
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


def contribution_values(contrib: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for panel, spec in PANELS.items():
        for region, items in [("CRE1", spec["left"]), ("CRE4", spec["right"])]:
            for label, gene in items:
                sub = contrib[(contrib["gene"].eq(gene)) & (contrib["region"].eq(region))]
                for tf in TF_STACK:
                    tf_sub = sub[sub["tf"].eq(tf)].sort_values(["seq_start", "seq_end", "site_index"])
                    if tf_sub.empty:
                        rows.append({
                            "panel": panel,
                            "base": spec["base"],
                            "region": region,
                            "label": label,
                            "gene": gene,
                            "tf": tf,
                            "site_key": f"{tf}:none",
                            "site_index": -1,
                            "seq_start": -1,
                            "seq_end": -1,
                            "activation_contribution": 0.0,
                        })
                        continue
                    for row in tf_sub.itertuples(index=False):
                        rows.append({
                            "panel": panel,
                            "base": spec["base"],
                            "region": region,
                            "label": label,
                            "gene": gene,
                            "tf": tf,
                            "site_key": f"{tf}:{int(row.seq_start)}-{int(row.seq_end)}:{int(row.site_index)}",
                            "site_index": int(row.site_index),
                            "seq_start": int(row.seq_start),
                            "seq_end": int(row.seq_end),
                            "activation_contribution": float(row.activation_contribution),
                        })
    return pd.DataFrame(rows)


def draw_panel(body: list[str], values: pd.DataFrame, panel: str, x: float, y: float, w: float, h: float) -> None:
    spec = PANELS[panel]
    groups = [("CRE1", *item) for item in spec["left"]] + [("CRE4", *item) for item in spec["right"]]
    axis_x = x + 42
    axis_top = y + 28
    axis_bottom = y + h - 35
    axis_right = x + w - 8
    ymax = 1.85
    header_h = 18
    bar_w = 28

    body.append(text(x - 25, y - 5, panel, 19, "700"))
    body.append(text(axis_x + 8, y + 34, f"→{spec['base']}", 20, "700"))

    cre1_l = axis_x + 24
    cre1_r = axis_x + (axis_right - axis_x) * 0.49
    cre4_l = axis_x + (axis_right - axis_x) * 0.53
    cre4_r = axis_right
    body.append(rect(cre1_l, y, cre1_r - cre1_l, header_h, "#D0D0D0", "none"))
    body.append(rect(cre4_l, y, cre4_r - cre4_l, header_h, "#D0D0D0", "none"))
    body.append(text((cre1_l + cre1_r) / 2, y + 14, "CRE1", 13, anchor="middle"))
    body.append(text((cre4_l + cre4_r) / 2, y + 14, "CRE4", 13, anchor="middle"))

    body.append(line(axis_x, axis_top, axis_x, axis_bottom, "#222222", 1.1))
    body.append(line(axis_x, axis_bottom, axis_right, axis_bottom, "#222222", 1.1))
    for tick in [0, 0.5, 1.0, 1.5]:
        ty = scale(tick, 0, ymax, axis_bottom, axis_top)
        body.append(line(axis_x - 5, ty, axis_x, ty, "#222222", 1.0))
        body.append(text(axis_x - 8, ty + 4, f"{tick:.1f}", 12, anchor="end"))
    body.append(text(x - 20, (axis_top + axis_bottom) / 2, "ΔΔA", 16, anchor="middle", extra=f'transform="rotate(-90 {x - 20:.2f} {(axis_top + axis_bottom) / 2:.2f})"'))

    x_positions = []
    for i in range(8):
        if i < 4:
            bx = scale(i, 0, 3, cre1_l + 26, cre1_r - 26)
        else:
            bx = scale(i - 4, 0, 3, cre4_l + 26, cre4_r - 26)
        x_positions.append(bx)

    boundaries: list[list[float]] = []
    for i, (region, label, gene) in enumerate(groups):
        bx = x_positions[i]
        bottom = axis_bottom
        bounds = [axis_bottom]
        sub = values[(values["panel"].eq(panel)) & (values["region"].eq(region)) & (values["gene"].eq(gene))]
        for tf in TF_STACK:
            tf_sub = sub[sub["tf"].eq(tf)].sort_values(["seq_start", "seq_end", "site_index"])
            for row in tf_sub.itertuples(index=False):
                val = float(row.activation_contribution)
                if val <= 0:
                    continue
                bh = scale(val, 0, ymax, 0, axis_bottom - axis_top)
                if bh > 0.25:
                    body.append(rect(bx - bar_w / 2, bottom - bh, bar_w, bh, TF_COLORS[tf], "#222222", 0.55, 0.95))
                    body.append(line(bx - bar_w / 2, bottom - bh, bx + bar_w / 2, bottom - bh, "#222222", 0.45))
                bottom -= bh
                bounds.append(bottom)
        boundaries.append(bounds)
        body.append(text(bx, y + h - 10, label, 11, anchor="middle"))

    # Connect cumulative site-level segment boundaries between adjacent variants.
    for start, end in [(0, 3), (4, 7)]:
        for i in range(start, end):
            n = min(len(boundaries[i]), len(boundaries[i + 1]))
            for j in range(1, n):
                y1 = boundaries[i][j]
                y2 = boundaries[i + 1][j]
                body.append(line(x_positions[i] + bar_w / 2, y1, x_positions[i + 1] - bar_w / 2, y2, "#333333", 0.75))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    genes = sorted({gene for spec in PANELS.values() for _, gene in spec["left"] + spec["right"]})
    coefs = fig5.coefficients(XML)
    contrib = fig5.contribution_table(XML, genes, coefs)
    contrib.to_csv(OUT / "fig5_fk_site_level_contributions.csv", index=False)
    values = contribution_values(contrib)
    values.to_csv(OUT / "fig5_fk_stacked_bar_values.csv", index=False)

    width, height = 880, 235
    body: list[str] = []
    draw_panel(body, values, "F", 36, 20, 385, 190)
    draw_panel(body, values, "K", 477, 20, 385, 190)

    write_svg(OUT / "fig5_fk_paper_exact.svg", width, height, body)
    write_svg(OUT / "fig5_fk_paper_exact_4x.svg", width, height, body, scale_factor=4)
    print(OUT / "fig5_fk_paper_exact.svg")


if __name__ == "__main__":
    main()
