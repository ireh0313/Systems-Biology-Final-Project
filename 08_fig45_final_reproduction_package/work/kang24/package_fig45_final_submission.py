#!/usr/bin/env python3
"""Assemble final submission-ready Fig. 4/Fig. 5 plates and notes."""

from __future__ import annotations

import csv
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "outputs" / "fig45_competition_off_repeats"
FINAL = ROOT / "outputs" / "fig45_final_submission"
SEED = "4tf_seed0"


def copy_panel(src: Path, dest_name: str) -> str:
    FINAL.mkdir(parents=True, exist_ok=True)
    dest = FINAL / dest_name
    shutil.copyfile(src, dest)
    return dest.name


def read_rows() -> list[dict[str, str]]:
    with (RUN / "fig45_parallel_repeats_summary.csv").open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def fnum(value: str) -> float:
    return float(value) if value else float("nan")


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def median(values: list[float]) -> float:
    values = sorted(values)
    n = len(values)
    mid = n // 2
    return values[mid] if n % 2 else (values[mid - 1] + values[mid]) / 2


def image(file_name: str, x: int, y: int, w: int, h: int) -> str:
    return (
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="#fff" '
        f'stroke="#d7d7d7" stroke-width="1"/>'
        f'<image href="{file_name}" x="{x + 8}" y="{y + 8}" '
        f'width="{w - 16}" height="{h - 16}" preserveAspectRatio="xMidYMid meet"/>'
    )


def label(text: str, x: int, y: int, size: int = 28, weight: int = 700) -> str:
    return (
        f'<text x="{x}" y="{y}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" fill="#222">{text}</text>'
    )


def caption(text: str, x: int, y: int, size: int = 20) -> str:
    return (
        f'<text x="{x}" y="{y}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" fill="#444">{text}</text>'
    )


def write_svg(path: Path, width: int, height: int, body: list[str]) -> None:
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        *body,
        "</svg>",
    ]
    path.write_text("\n".join(svg) + "\n", encoding="utf-8")


def main() -> None:
    rows = read_rows()
    seeds = [r for r in rows if r["kind"] == "seed"]
    scrambles = [r for r in rows if r["kind"] == "scramble"]
    seed0 = next(r for r in seeds if r["id"] == "0")

    seed_rs = [fnum(r["singlehit_r"]) for r in seeds]
    fig4j_rs = [fnum(r["fig4j_multihit_r"]) for r in seeds]
    fig5_rs = [fnum(r["fig5_at_r"]) for r in seeds]
    scramble_rs = [fnum(r["singlehit_r_against_true_labels"]) for r in scrambles]

    fig4_panels = {
        "profile": copy_panel(
            RUN / SEED / "kang24_4tf_seed0_singlehit_delta_profile.svg",
            "fig4_A_seed0_singlehit_profile.svg",
        ),
        "singlehit": copy_panel(
            RUN / SEED / "kang24_4tf_seed0_singlehit_delta_scatter.svg",
            "fig4_B_seed0_singlehit_scatter.svg",
        ),
        "multihit": copy_panel(
            RUN / SEED / "fig4j" / "kang24_fig4j_4tf_multihit_scatter.svg",
            "fig4_C_seed0_multihit_scatter.svg",
        ),
        "scramble": copy_panel(
            RUN / "fig45_seed_vs_scramble_correlation.svg",
            "fig4_D_seed_vs_scramble_correlation.svg",
        ),
    }

    fig5_panels = {
        "at": copy_panel(
            RUN / SEED / "fig5" / "kang24_fig5_AT_substitution_profile.svg",
            "fig5_A_AT_substitution_profile.svg",
        ),
        "occupancy": copy_panel(
            RUN / SEED / "fig5" / "kang24_fig5_TFBS_occupancy_maps.svg",
            "fig5_B_TFBS_occupancy_maps.svg",
        ),
        "coefficients": copy_panel(
            RUN / SEED / "fig5" / "kang24_fig5_activation_coefficients.svg",
            "fig5_C_activation_coefficients.svg",
        ),
        "contrib": copy_panel(
            RUN / SEED / "fig5" / "kang24_fig5_CRE1_CRE4_contribution_bars.svg",
            "fig5_D_CRE1_CRE4_contribution_bars.svg",
        ),
    }

    common = (
        "4TF_self full fitting; Competition=false; SelfCompetition=true; "
        "representative seed=0; 5 seeds and 10 scrambled-label controls"
    )

    fig4_body = [
        label("Figure 4 reproduction package", 48, 56, 34),
        caption(common, 48, 88, 18),
        label("A", 48, 140),
        caption("Single-hit profile, seed 0", 86, 140),
        image(fig4_panels["profile"], 48, 160, 1480, 620),
        label("B", 48, 850),
        caption(f"Single-hit fit, r={fnum(seed0['singlehit_r']):.3f}", 86, 850),
        image(fig4_panels["singlehit"], 48, 870, 700, 620),
        label("C", 828, 850),
        caption(f"Multi-hit validation, r={fnum(seed0['fig4j_multihit_r']):.3f}", 866, 850),
        image(fig4_panels["multihit"], 828, 870, 700, 620),
        label("D", 48, 1580),
        caption(
            f"Real labels vs scramble controls: real mean r={mean(seed_rs):.3f}, scramble mean r={mean(scramble_rs):.3f}",
            86,
            1580,
        ),
        image(fig4_panels["scramble"], 48, 1600, 1480, 760),
        caption(
            f"Across real seeds: single-hit mean={mean(seed_rs):.3f}, median={median(seed_rs):.3f}; "
            f"Fig. 4J multi-hit mean={mean(fig4j_rs):.3f}.",
            48,
            2430,
            18,
        ),
    ]
    write_svg(FINAL / "final_fig4_competition_off.svg", 1580, 2480, fig4_body)

    fig5_body = [
        label("Figure 5 reproduction package", 48, 56, 34),
        caption(common, 48, 88, 18),
        label("A", 48, 140),
        caption(f"A/T substitution activity profile, seed 0, r={fnum(seed0['fig5_at_r']):.3f}", 86, 140),
        image(fig5_panels["at"], 48, 160, 1480, 650),
        label("B", 48, 890),
        caption("Predicted TFBS occupancy maps", 86, 890),
        image(fig5_panels["occupancy"], 48, 920, 1480, 760),
        label("C", 48, 1760),
        caption("Activation coefficients", 86, 1760),
        image(fig5_panels["coefficients"], 48, 1780, 700, 520),
        label("D", 828, 1760),
        caption("CRE1/CRE4 contribution proxy", 866, 1760),
        image(fig5_panels["contrib"], 828, 1780, 700, 620),
        caption(
            f"Across real seeds: Fig. 5 A/T profile mean r={mean(fig5_rs):.3f}, "
            f"median={median(fig5_rs):.3f}, range={min(fig5_rs):.3f}-{max(fig5_rs):.3f}.",
            48,
            2470,
            18,
        ),
    ]
    write_svg(FINAL / "final_fig5_competition_off.svg", 1580, 2520, fig5_body)

    (FINAL / "fig45_final_metrics.csv").write_text(
        "\n".join(
            [
                "metric,value",
                f"representative_seed,{SEED}",
                f"seed0_singlehit_r,{fnum(seed0['singlehit_r']):.6f}",
                f"seed0_fig4j_multihit_r,{fnum(seed0['fig4j_multihit_r']):.6f}",
                f"seed0_fig5_at_r,{fnum(seed0['fig5_at_r']):.6f}",
                f"real_seed_singlehit_mean,{mean(seed_rs):.6f}",
                f"real_seed_singlehit_median,{median(seed_rs):.6f}",
                f"real_seed_fig4j_mean,{mean(fig4j_rs):.6f}",
                f"real_seed_fig5_mean,{mean(fig5_rs):.6f}",
                f"scramble_true_label_r_mean,{mean(scramble_rs):.6f}",
                f"scramble_true_label_r_median,{median(scramble_rs):.6f}",
                f"scramble_true_label_r_min,{min(scramble_rs):.6f}",
                f"scramble_true_label_r_max,{max(scramble_rs):.6f}",
                f"scramble_repeats,{len(scramble_rs)}",
                f"real_seed_repeats,{len(seed_rs)}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(FINAL)


if __name__ == "__main__":
    main()
