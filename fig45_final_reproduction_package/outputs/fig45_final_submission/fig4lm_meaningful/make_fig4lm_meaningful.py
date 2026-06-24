#!/usr/bin/env python3
"""Create a meaning-focused Figure 4L/M replacement from available results."""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "outputs" / "fig45_competition_off_repeats"
SUMMARY = RUN / "fig45_parallel_repeats_summary.csv"
SCR_MULTI = ROOT / "outputs" / "fig45_final_submission" / "fig4lm_paper_exact" / "scramble_multihit" / "scramble_multihit_correlations.csv"
OUT = ROOT / "outputs" / "fig45_final_submission" / "fig4lm_meaningful"

RED = "#F8766D"
CYAN = "#00BFC4"
DARK = "#222222"
GRID = "#E7E7E7"


def esc(s: object) -> str:
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def scale(v: float, a: float, b: float, c: float, d: float) -> float:
    if a == b:
        return (c + d) / 2
    return c + (v - a) / (b - a) * (d - c)


def text(x: float, y: float, s: str, size: int = 16, weight: str = "400", anchor: str = "start", fill: str = DARK) -> str:
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" fill="{fill}">{esc(s)}</text>'
    )


def vtext(x: float, y: float, s: str, size: int = 19) -> str:
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-family="Arial, Helvetica, sans-serif" font-size="{size}" '
        f'text-anchor="middle" fill="{DARK}" transform="rotate(-90 {x:.2f} {y:.2f})">{esc(s)}</text>'
    )


def line(x1: float, y1: float, x2: float, y2: float, color: str = DARK, width: float = 1.0, dash: str | None = None) -> str:
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{color}" stroke-width="{width}"{dash_attr}/>'


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


def mean_se(values: list[float]) -> tuple[float, float]:
    series = pd.Series(values, dtype=float)
    return float(series.mean()), float(series.std(ddof=1) / math.sqrt(len(series)))


def jitter_positions(n: int, center: float, spread: float) -> list[float]:
    if n == 1:
        return [center]
    return [center + (i - (n - 1) / 2) * spread / max(1, n - 1) for i in range(n)]


def panel(letter: str, subtitle: str, real_values: list[float], control_values: list[float], x: float, y: float, w: float, h: float, show_ylabel: bool) -> list[str]:
    ylo, yhi = -0.7, 0.9
    body: list[str] = [text(x - 34, y - 18, letter, 20, "700")]
    body.append(rect(x, y, w, h, "#FFFFFF"))
    for tick in [-0.6, -0.3, 0.0, 0.3, 0.6, 0.9]:
        ty = scale(tick, ylo, yhi, y + h, y)
        body.append(line(x, ty, x + w, ty, GRID, 1.0))
        if show_ylabel:
            body.append(text(x - 9, ty + 5, f"{tick:.1f}", 14, anchor="end", fill="#4c4c4c"))
    body.append(line(x, scale(0, ylo, yhi, y + h, y), x + w, scale(0, ylo, yhi, y + h, y), "#999999", 1.2, "4 4"))
    body.append(line(x, y + h, x + w, y + h, DARK, 1.2))
    body.append(line(x, y, x, y + h, DARK, 1.2))

    centers = {"real": x + w * 0.34, "control": x + w * 0.70}
    for label, values, color in [("real", real_values, RED), ("control", control_values, CYAN)]:
        center = centers[label]
        xs = jitter_positions(len(values), center, 36)
        for px, value in zip(xs, values):
            py = scale(value, ylo, yhi, y + h, y)
            body.append(f'<circle cx="{px:.2f}" cy="{py:.2f}" r="4.2" fill="{color}" opacity="0.82"/>')
        mean, se = mean_se(values)
        mean_y = scale(mean, ylo, yhi, y + h, y)
        se_low = scale(mean - se, ylo, yhi, y + h, y)
        se_high = scale(mean + se, ylo, yhi, y + h, y)
        body.append(line(center, se_low, center, se_high, color, 2.0))
        body.append(line(center - 15, se_low, center + 15, se_low, color, 2.0))
        body.append(line(center - 15, se_high, center + 15, se_high, color, 2.0))
        body.append(line(center - 24, mean_y, center + 24, mean_y, color, 4.0))
        label_y = mean_y + 18 if mean > 0.7 else mean_y - 12
        body.append(text(center, label_y, f"{mean:.2f}", 12, "700", anchor="middle", fill=color))

    body.append(text(centers["real"], y + h + 24, "4TF", 15, anchor="middle"))
    body.append(text(centers["control"], y + h + 24, "scramble", 15, anchor="middle"))
    body.append(text(x + w / 2, y - 4, subtitle, 13, anchor="middle", fill="#333333"))
    if show_ylabel:
        body.append(vtext(x - 58, y + h / 2, "Pearson's R", 19))
    return body


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = pd.read_csv(SUMMARY)
    scr_multi = pd.read_csv(SCR_MULTI)

    seeds = summary[summary["kind"].eq("seed")].copy()
    scrambles = summary[summary["kind"].eq("scramble")].copy()
    l_real = pd.to_numeric(seeds["singlehit_r"], errors="coerce").dropna().tolist()
    l_ctrl = pd.to_numeric(scrambles["singlehit_r_against_true_labels"], errors="coerce").dropna().tolist()
    m_real = pd.to_numeric(seeds["fig4j_multihit_r"], errors="coerce").dropna().tolist()
    m_ctrl = pd.to_numeric(scr_multi["multihit_r_against_true_labels"], errors="coerce").dropna().tolist()

    width, height = 780, 440
    body: list[str] = []
    body.append(text(width / 2, 28, "4TF model performance against scrambled-label controls", 20, "400", anchor="middle"))
    body += panel("L", "single-hit training set", l_real, l_ctrl, 105, 72, 260, 265, True)
    body += panel("M", "multi-hit validation set", m_real, m_ctrl, 455, 72, 260, 265, False)
    body.append(text(390, 407, "Each point is one independent fit/control; thick line = mean, whiskers = SE.", 13, anchor="middle", fill="#555555"))

    svg = OUT / "fig4lm_meaningful_4tf_vs_scramble.svg"
    svg_4x = OUT / "fig4lm_meaningful_4tf_vs_scramble_4x.svg"
    write_svg(svg, width, height, body)
    write_svg(svg_4x, width, height, body, scale_factor=4)

    rows = []
    for panel_id, series, values in [
        ("L", "4TF real seeds", l_real),
        ("L", "scrambled-label controls", l_ctrl),
        ("M", "4TF real seeds", m_real),
        ("M", "scrambled-label controls", m_ctrl),
    ]:
        mean, se = mean_se(values)
        rows.append({"panel": panel_id, "series": series, "mean": mean, "se": se, "n": len(values), "values": ";".join(f"{v:.6f}" for v in values)})
    pd.DataFrame(rows).to_csv(OUT / "fig4lm_meaningful_4tf_vs_scramble_metrics.csv", index=False)
    print(svg)


if __name__ == "__main__":
    main()
