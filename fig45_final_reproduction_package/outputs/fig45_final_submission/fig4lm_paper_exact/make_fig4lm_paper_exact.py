#!/usr/bin/env python3
"""Create paper-style Figure 4L/M panels using only reproduced values."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
KANG = ROOT / "work" / "kang24"
sys.path.insert(0, str(KANG))
import make_fig4j_multihit as fig4j  # noqa: E402


RUN = ROOT / "outputs" / "fig45_competition_off_repeats"
SUMMARY = RUN / "fig45_parallel_repeats_summary.csv"
OUT = ROOT / "outputs" / "fig45_final_submission" / "fig4lm_paper_exact"

RED = "#F8766D"
CYAN = "#00BFC4"
AXIS = "#222222"
GRID = "#E6E6E6"


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


def vtext(x: float, y: float, s: str, size: int = 20) -> str:
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-family="Arial, Helvetica, sans-serif" font-size="{size}" '
        f'text-anchor="middle" fill="#111111" transform="rotate(-90 {x:.2f} {y:.2f})">{esc(s)}</text>'
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


def mean_se(values: pd.Series) -> tuple[float, float, int]:
    vals = pd.to_numeric(values, errors="coerce").dropna()
    se = float(vals.std(ddof=1) / math.sqrt(len(vals))) if len(vals) > 1 else 0.0
    return float(vals.mean()), se, int(len(vals))


def scramble_multihit_summary(summary: pd.DataFrame) -> pd.DataFrame:
    out_dir = OUT / "scramble_multihit"
    out_dir.mkdir(parents=True, exist_ok=True)
    cache = out_dir / "scramble_multihit_correlations.csv"
    if cache.exists():
        return pd.read_csv(cache)

    seqs = fig4j.read_multihit("multi-hit")
    rows = []
    for row in summary[summary["kind"].eq("scramble")].itertuples(index=False):
        xml = Path(row.xml)
        work_xml = out_dir / f"{xml.stem}_multihit.xml"
        fig4j.build_multihit_xml(xml, seqs, work_xml)
        predicted = fig4j.run_rate(work_xml)
        merged = fig4j.add_delta(seqs.merge(predicted, on="xml_name", how="inner"))
        merged.to_csv(out_dir / f"{xml.stem}_multihit_observed_predicted.csv", index=False)
        rows.append(
            {
                "id": row.id,
                "kind": "scramble",
                "multihit_r_against_true_labels": merged["observed_delta"].corr(merged["predicted_delta"]),
                "n": len(merged),
                "xml": str(xml),
            }
        )
    result = pd.DataFrame(rows).sort_values("id")
    result.to_csv(cache, index=False)
    return result


def y_map(value: float, x: float, y: float, top_h: float, bot_h: float, gap: float) -> float:
    top_y = y
    bot_y = y + top_h + gap
    top_range = (0.68, 0.86)
    bot_range = (0.00, 0.55)
    if value >= top_range[0]:
        return scale(value, top_range[0], top_range[1], top_y + top_h, top_y)
    return scale(max(value, bot_range[0]), bot_range[0], bot_range[1], bot_y + bot_h, bot_y)


def error_point(body: list[str], px: float, py: float, mean: float, se: float, x: float, y: float, top_h: float, bot_h: float, gap: float, color: str) -> None:
    y1 = y_map(mean - se, x, y, top_h, bot_h, gap)
    y2 = y_map(mean + se, x, y, top_h, bot_h, gap)
    if mean >= 0:
        body.append(line(px, y1, px, y2, color, 1.8))
        body.append(line(px - 7, y1, px + 7, y1, color, 1.8))
        body.append(line(px - 7, y2, px + 7, y2, color, 1.8))
        body.append(f'<circle cx="{px:.2f}" cy="{py:.2f}" r="4.6" fill="{color}" stroke="{color}" stroke-width="1"/>')
    else:
        body.append(f'<path d="M {px-5:.2f} {py-2:.2f} L {px+5:.2f} {py-2:.2f} L {px:.2f} {py+8:.2f} Z" fill="{color}" opacity="0.95"/>')
        body.append(line(px, py - 18, px, py - 4, color, 1.5))


def panel(letter: str, real: tuple[float, float, int], control: tuple[float, float, int], x: float, y: float, w: float, h: float, show_ylabel: bool, show_legend: bool) -> list[str]:
    body: list[str] = [text(x - 42, y - 24, letter, 20, "700")]
    top_h = 70
    gap = 18
    bot_h = h - top_h - gap
    plot_x = x
    plot_w = w
    top_y = y
    bot_y = y + top_h + gap

    body.append(rect(plot_x, top_y, plot_w, top_h, "#ffffff"))
    body.append(rect(plot_x, bot_y, plot_w, bot_h, "#ffffff"))
    xmin, xmax = 1.75, 5.25
    for tick in [2, 3, 4, 5]:
        tx = scale(tick, xmin, xmax, plot_x, plot_x + plot_w)
        body.append(line(tx, top_y, tx, top_y + top_h, GRID, 0.9))
        body.append(line(tx, bot_y, tx, bot_y + bot_h, GRID, 0.9))
        body.append(line(tx, bot_y + bot_h, tx, bot_y + bot_h + 6, AXIS, 1.1))
        body.append(text(tx, bot_y + bot_h + 28, str(tick), 18, anchor="middle", fill="#4c4c4c"))
    for tick in [0.7, 0.8]:
        ty = y_map(tick, x, y, top_h, bot_h, gap)
        body.append(line(plot_x, ty, plot_x + plot_w, ty, GRID, 1.0))
        if show_ylabel:
            body.append(text(plot_x - 10, ty + 6, f"{tick:.1f}", 18, anchor="end", fill="#4c4c4c"))
    for tick in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]:
        ty = y_map(tick, x, y, top_h, bot_h, gap)
        body.append(line(plot_x, ty, plot_x + plot_w, ty, GRID, 1.0))
        if show_ylabel:
            body.append(text(plot_x - 10, ty + 6, f"{tick:.1f}", 18, anchor="end", fill="#4c4c4c"))

    for yy, hh in [(top_y, top_h), (bot_y, bot_h)]:
        body.append(line(plot_x, yy, plot_x, yy + hh, AXIS, 1.2))
        body.append(line(plot_x, yy + hh, plot_x + plot_w, yy + hh, AXIS, 1.2))
    body.append(line(plot_x - 8, top_y + top_h - 3, plot_x + 8, top_y + top_h - 15, AXIS, 1.2))
    body.append(line(plot_x - 8, bot_y + 15, plot_x + 8, bot_y + 3, AXIS, 1.2))

    x4 = scale(4, xmin, xmax, plot_x, plot_x + plot_w)
    real_mean, real_se, _ = real
    ctrl_mean, ctrl_se, _ = control
    error_point(body, x4, y_map(real_mean, x, y, top_h, bot_h, gap), real_mean, real_se, x, y, top_h, bot_h, gap, RED)
    error_point(body, x4, y_map(ctrl_mean, x, y, top_h, bot_h, gap), ctrl_mean, ctrl_se, x, y, top_h, bot_h, gap, CYAN)
    if ctrl_mean < 0:
        body.append(text(x4, bot_y + bot_h - 4, f"{ctrl_mean:.2f}", 10, anchor="middle", fill=CYAN))

    body.append(text(plot_x + plot_w / 2, bot_y + bot_h + 54, "Number of TFs", 22, anchor="middle"))
    if show_ylabel:
        body.append(vtext(plot_x - 68, y + h / 2, "Pearson's R", 21))
    if show_legend:
        lx = plot_x + 18
        ly = bot_y + bot_h - 54
        body.append(line(lx, ly, lx + 18, ly, RED, 2.2))
        body.append(f'<circle cx="{lx + 9:.2f}" cy="{ly:.2f}" r="4.2" fill="{RED}"/>')
        body.append(text(lx + 24, ly + 6, "4TF_self real seeds", 14))
        body.append(line(lx, ly + 24, lx + 18, ly + 24, CYAN, 2.2))
        body.append(f'<circle cx="{lx + 9:.2f}" cy="{ly + 24:.2f}" r="4.2" fill="{CYAN}"/>')
        body.append(text(lx + 24, ly + 30, "scrambled labels", 14))
    return body


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = pd.read_csv(SUMMARY)
    scr_multi = scramble_multihit_summary(summary)

    seeds = summary[summary["kind"].eq("seed")]
    scrambles = summary[summary["kind"].eq("scramble")]
    l_real = mean_se(seeds["singlehit_r"])
    l_control = mean_se(scrambles["singlehit_r_against_true_labels"])
    m_real = mean_se(seeds["fig4j_multihit_r"])
    m_control = mean_se(scr_multi["multihit_r_against_true_labels"])

    width, height = 610, 390
    body: list[str] = []
    body += panel("L", l_real, l_control, 102, 38, 190, 270, True, True)
    body += panel("M", m_real, m_control, 352, 38, 190, 270, False, False)

    svg = OUT / "fig4lm_our_values_paper_exact.svg"
    svg_4x = OUT / "fig4lm_our_values_paper_exact_4x.svg"
    write_svg(svg, width, height, body)
    write_svg(svg_4x, width, height, body, scale_factor=4)

    pd.DataFrame(
        [
            {"panel": "L", "series": "4TF_self real seeds", "mean": l_real[0], "se": l_real[1], "n": l_real[2]},
            {"panel": "L", "series": "scrambled-label controls", "mean": l_control[0], "se": l_control[1], "n": l_control[2]},
            {"panel": "M", "series": "4TF_self real seeds", "mean": m_real[0], "se": m_real[1], "n": m_real[2]},
            {"panel": "M", "series": "scrambled-label controls", "mean": m_control[0], "se": m_control[1], "n": m_control[2]},
        ]
    ).to_csv(OUT / "fig4lm_our_values_paper_exact_metrics.csv", index=False)
    print(svg)


if __name__ == "__main__":
    main()
