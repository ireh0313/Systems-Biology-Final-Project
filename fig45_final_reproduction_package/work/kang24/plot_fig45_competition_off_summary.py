#!/usr/bin/env python3
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SUMMARY = ROOT / "outputs/fig45_competition_off_repeats/fig45_parallel_repeats_summary.csv"
OUT = ROOT / "outputs/fig45_competition_off_repeats/fig45_seed_vs_scramble_correlation.svg"


def mean(values):
    return sum(values) / len(values)


def median(values):
    ordered = sorted(values)
    n = len(ordered)
    mid = n // 2
    if n % 2:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2


with SUMMARY.open(newline="") as f:
    rows = list(csv.DictReader(f))

seed = [float(r["singlehit_r"]) for r in rows if r["kind"] == "seed"]
scramble = [
    float(r["singlehit_r_against_true_labels"])
    for r in rows
    if r["kind"] == "scramble"
]

w, h = 920, 520
margin = 72
plot_w = w - margin * 2
plot_h = h - 140
y_min, y_max = -0.7, 0.9


def y(v):
    return margin + (y_max - v) / (y_max - y_min) * plot_h


def x(group, i, n):
    center = margin + plot_w * (0.30 if group == "seed" else 0.70)
    if n == 1:
        return center
    return center + (i - (n - 1) / 2) * 20


def line(x1, y1, x2, y2, stroke="#222", width=1, dash=None):
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    return (
        f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
        f'stroke="{stroke}" stroke-width="{width}"{dash_attr}/>'
    )


parts = [
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
    '<rect width="100%" height="100%" fill="#ffffff"/>',
    '<text x="72" y="38" font-family="Arial" font-size="22" font-weight="700">Corrected Fig. 4/5 4TF fit: real labels vs scrambled labels</text>',
    '<text x="72" y="64" font-family="Arial" font-size="13" fill="#444">Competition=false, SelfCompetition=true; 5 real-seed fits and 10 scrambled-label controls</text>',
]

for tick in [-0.6, -0.3, 0.0, 0.3, 0.6, 0.9]:
    yy = y(tick)
    parts.append(line(margin, yy, w - margin, yy, "#d8d8d8", 1))
    parts.append(
        f'<text x="42" y="{yy + 4:.1f}" font-family="Arial" font-size="12" fill="#555">{tick:.1f}</text>'
    )

parts.append(line(margin, y(0), w - margin, y(0), "#777", 1.4, "5,5"))
parts.append(
    f'<text x="30" y="{margin + plot_h / 2:.1f}" transform="rotate(-90 30 {margin + plot_h / 2:.1f})" font-family="Arial" font-size="13" fill="#333">Pearson r vs observed single-hit activity</text>'
)

seed_color = "#1f77b4"
scramble_color = "#c44e52"
for i, value in enumerate(seed):
    xx, yy = x("seed", i, len(seed)), y(value)
    parts.append(f'<circle cx="{xx:.1f}" cy="{yy:.1f}" r="7" fill="{seed_color}" opacity="0.92"/>')
for i, value in enumerate(scramble):
    xx, yy = x("scramble", i, len(scramble)), y(value)
    parts.append(f'<circle cx="{xx:.1f}" cy="{yy:.1f}" r="7" fill="{scramble_color}" opacity="0.88"/>')

seed_mean = mean(seed)
scramble_mean = mean(scramble)
parts.append(line(x("seed", 0, 1) - 80, y(seed_mean), x("seed", 0, 1) + 80, y(seed_mean), seed_color, 4))
parts.append(line(x("scramble", 0, 1) - 110, y(scramble_mean), x("scramble", 0, 1) + 110, y(scramble_mean), scramble_color, 4))

for group, values, color, label in [
    ("seed", seed, seed_color, "Real labels"),
    ("scramble", scramble, scramble_color, "Scrambled labels"),
]:
    cx = x(group, 0, 1)
    parts.append(
        f'<text x="{cx:.1f}" y="{h - 78}" font-family="Arial" font-size="16" text-anchor="middle" font-weight="700" fill="#222">{label}</text>'
    )
    parts.append(
        f'<text x="{cx:.1f}" y="{h - 54}" font-family="Arial" font-size="13" text-anchor="middle" fill="#444">n={len(values)}, mean={mean(values):.3f}, median={median(values):.3f}</text>'
    )
    parts.append(
        f'<text x="{cx:.1f}" y="{h - 34}" font-family="Arial" font-size="13" text-anchor="middle" fill="#444">range {min(values):.3f} to {max(values):.3f}</text>'
    )

parts.append("</svg>")
OUT.write_text("\n".join(parts) + "\n")
print(OUT)
