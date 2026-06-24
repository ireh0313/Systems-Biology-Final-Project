from pathlib import Path
import csv
import math

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D
import numpy as np


ROOT = Path.home() / "cre_repro"
PWM_CSV = ROOT / "work" / "csv" / "pwms.csv"
OUTPUT = ROOT / "results" / "fig2" / "paper_method"
TF_ORDER = ["CREB1", "CREB3", "CREB5", "CREM", "ATF1", "ATF4", "ATF7"]
BASES = ["A", "C", "G", "T"]
COLORS = {
    "A": "#009B72",
    "C": "#28689B",
    "G": "#F2A719",
    "T": "#D9233F",
}
FONT = FontProperties(family="DejaVu Sans", weight="bold")


def read_log_odds(path):
    matrices = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.reader(handle))

    index = 0
    while index < len(rows):
        row = rows[index]
        if row and row[0] in TF_ORDER and len(row) > 2:
            tf_name = row[0]
            block = rows[index : index + 4]
            values = {}
            for block_row in block:
                base = block_row[1]
                values[base] = [float(value) for value in block_row[2:] if value != ""]
            matrices[tf_name] = np.array([values[base] for base in BASES]).T
            index += 4
        else:
            index += 1
    return matrices


def log_odds_to_bits(scores):
    # Table S2 values are base-2 log-odds scores. With equal background
    # frequencies, probabilities are proportional to 0.25 * 2**score.
    weights = 0.25 * np.power(2.0, scores)
    probabilities = weights / weights.sum(axis=1, keepdims=True)
    entropy = -np.sum(
        np.where(probabilities > 0, probabilities * np.log2(probabilities), 0),
        axis=1,
    )
    information = 2.0 - entropy
    return probabilities * information[:, None]


def draw_letter(ax, letter, x, y, height):
    if height <= 0.002:
        return

    path = TextPath((0, 0), letter, size=1, prop=FONT)
    bounds = path.get_extents()
    target_width = 0.92
    sx = target_width / bounds.width
    sy = height / bounds.height
    transform = (
        Affine2D()
        .translate(-bounds.x0, -bounds.y0)
        .scale(sx, sy)
        .translate(x - target_width / 2, y)
    )
    ax.add_patch(
        PathPatch(
            path,
            transform=transform + ax.transData,
            facecolor=COLORS[letter],
            edgecolor="none",
        )
    )


def draw_logo(ax, heights, name, title_x):
    for position, row in enumerate(heights, start=1):
        bottom = 0.0
        for value, base in sorted(zip(row, BASES)):
            draw_letter(ax, base, position, bottom, value)
            bottom += value

    length = heights.shape[0]
    ax.set_xlim(0.45, length + 0.55)
    ax.set_ylim(0, 2.05)
    ax.set_xticks([])
    ax.set_yticks([])

    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_linewidth(1.35)
    ax.spines["bottom"].set_linewidth(1.35)
    ax.tick_params(length=0)
    ax.text(
        title_x,
        2.16,
        name,
        fontsize=20,
        fontweight="normal",
        ha="left",
        va="bottom",
        clip_on=False,
    )


def main():
    matrices = read_log_odds(PWM_CSV)
    missing = [name for name in TF_ORDER if name not in matrices]
    if missing:
        raise RuntimeError(f"Missing PWM matrices: {missing}")

    # Every motif position uses the same physical width. Therefore x-axis
    # length varies with motif length while nucleotide letter size stays fixed.
    fig = plt.figure(figsize=(5.2, 15.4), facecolor="white")
    max_length = max(matrix.shape[0] for matrix in matrices.values())
    axis_left = 0.24
    max_axis_width = 0.68
    axis_height = 0.085
    top = 0.875
    vertical_step = 0.137
    shared_title_x = matrices["CREB1"].shape[0] / 4
    axes = []

    for index, name in enumerate(TF_ORDER):
        length = matrices[name].shape[0]
        width = max_axis_width * length / max_length
        bottom = top - index * vertical_step
        ax = fig.add_axes([axis_left, bottom, width, axis_height])
        draw_logo(ax, log_odds_to_bits(matrices[name]), name, shared_title_x)
        axes.append(ax)

    crem_axis = axes[TF_ORDER.index("CREM")]
    crem_box = crem_axis.get_position()
    fig.text(
        crem_box.x0 - 0.105,
        crem_box.y1,
        "Bits",
        rotation=90,
        va="top",
        ha="center",
        fontsize=17,
    )

    OUTPUT.mkdir(parents=True, exist_ok=True)
    png = OUTPUT / "fig2d_7family_pwm_sequence_logos_bp_scaled.png"
    pdf = OUTPUT / "fig2d_7family_pwm_sequence_logos_bp_scaled.pdf"
    fig.savefig(png, dpi=500, facecolor="white")
    fig.savefig(pdf, facecolor="white")
    plt.close(fig)

    summary = OUTPUT / "fig2d_7family_pwm_sequence_logo_summary.csv"
    with summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["TF", "motif_length", "consensus"])
        for name in TF_ORDER:
            probabilities = 0.25 * np.power(2.0, matrices[name])
            probabilities /= probabilities.sum(axis=1, keepdims=True)
            consensus = "".join(BASES[index] for index in probabilities.argmax(axis=1))
            writer.writerow([name, len(probabilities), consensus])

    print(png)
    print(pdf)
    print(summary)


if __name__ == "__main__":
    main()
