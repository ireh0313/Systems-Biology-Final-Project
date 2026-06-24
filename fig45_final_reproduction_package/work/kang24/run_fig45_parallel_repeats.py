#!/usr/bin/env python3
"""Run Fig. 4/5 4TF seed repeats and scrambled-label controls in parallel."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import math
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
KANG = ROOT / "work" / "kang24"
TRANSC = ROOT / "work" / "repos" / "transcpp"
OUT = ROOT / "outputs" / "fig45_competition_off_repeats"
WORK = KANG / "fig45_competition_off_repeats"
PYTHON = Path(sys.executable)


def run_logged(cmd: list[str], cwd: Path, log: Path) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    with log.open("w", encoding="utf-8") as handle:
        handle.write(f"$ {' '.join(cmd)}\n")
        handle.flush()
        result = subprocess.run(cmd, cwd=cwd, text=True, stdout=handle, stderr=subprocess.STDOUT)
        elapsed = time.monotonic() - start
        handle.write(f"\nreal {elapsed:.2f}\n")
        handle.write(f"returncode {result.returncode}\n")
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd)


def run_capture(cmd: list[str], cwd: Path) -> str:
    return subprocess.run(cmd, cwd=cwd, check=True, text=True, capture_output=True).stdout


def xml_has_output(xml: Path) -> bool:
    return xml.exists() and "<Output>" in xml.read_text(errors="ignore")


def log_verified(log: Path) -> bool:
    return log.exists() and "output verified!" in log.read_text(errors="ignore")


def parse_log(log: Path) -> dict[str, float]:
    text = log.read_text(errors="ignore") if log.exists() else ""
    patterns = {
        "initial_score": r"The initial score is\s+([0-9.eE+-]+)",
        "after_initial_score": r"The score is\s+([0-9.eE+-]+)\s+after initial moves",
        "final_score": r"The final score is\s+([0-9.eE+-]+)",
        "total_steps": r"Total steps is\s+([0-9.eE+-]+)",
        "runtime_seconds": r"\nreal\s+([0-9.eE+-]+)",
    }
    out: dict[str, float] = {}
    for key, pattern in patterns.items():
        match = re.search(pattern, text)
        if match:
            out[key] = float(match.group(1))
    return out


def pearson(path: Path, x: str, y: str) -> float:
    if not path.exists():
        return math.nan
    df = pd.read_csv(path)
    return float(df[x].corr(df[y]))


def build_4tf_xml(seed: int, xml: Path, force: bool = False) -> None:
    if xml.exists() and not force:
        return
    xml.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(PYTHON),
        str(KANG / "build_kang24_xml.py"),
        "--model",
        "4TF_self",
        "--seed",
        str(seed),
        "--init-loop",
        "100000",
        "--trainable",
        "--competition",
        "false",
        "--self-competition",
        "true",
        "--out",
        str(xml),
    ]
    run_capture(cmd, ROOT)


def make_singlehit_figures(xml: Path, prefix: Path, title: str, force: bool = False) -> float:
    csv_path = prefix.with_name(prefix.name + "_observed_predicted.csv")
    if not csv_path.exists() or force:
        cmd = [
            str(PYTHON),
            str(KANG / "make_kang24_figures.py"),
            "--xml",
            str(xml),
            "--section",
            "Output",
            "--prefix",
            str(prefix),
            "--title",
            title,
            "--delta-activity",
        ]
        subprocess.run(cmd, cwd=ROOT, check=True)
    return pearson(csv_path, "observed_delta", "predicted_delta")


def make_fig4j(xml: Path, outdir: Path, force: bool = False) -> float:
    summary = outdir / "kang24_fig4j_summary.md"
    if not summary.exists() or force:
        cmd = [str(PYTHON), str(KANG / "make_fig4j_multihit.py"), "--model-xml", str(xml), "--outdir", str(outdir)]
        subprocess.run(cmd, cwd=ROOT, check=True)
    csv_path = outdir / "kang24_fig4j_4tf_multihit_observed_predicted.csv"
    return pearson(csv_path, "observed_delta", "predicted_delta")


def make_fig5(xml: Path, outdir: Path, force: bool = False) -> float:
    summary = outdir / "kang24_fig5_summary.md"
    if not summary.exists() or force:
        cmd = [str(PYTHON), str(KANG / "make_fig5_analysis.py"), "--xml", str(xml), "--outdir", str(outdir)]
        subprocess.run(cmd, cwd=ROOT, check=True)
    csv_path = outdir / "kang24_fig5_AT_substitution_table.csv"
    return pearson(csv_path, "observed_delta", "predicted_delta")


def fit_xml(xml: Path, log: Path, force: bool = False) -> None:
    if xml_has_output(xml) and log_verified(log) and not force:
        return
    run_logged([str(TRANSC / "transcpp"), str(xml.resolve())], TRANSC, log)


def scramble_xml(template: Path, scrambled: Path, force: bool = False) -> None:
    if scrambled.exists() and not force:
        return
    scrambled.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(TRANSC / "scramble"),
        "--permute=RateData",
        "--by",
        "col",
        str(template.resolve()),
        str(scrambled.resolve()),
    ]
    run_capture(cmd, TRANSC)


def run_seed(seed: int, force: bool = False) -> dict[str, object]:
    slug = f"4tf_seed{seed}"
    xml = WORK / f"kang24_{slug}.xml"
    log = WORK / "logs" / f"{slug}.log"
    build_4tf_xml(seed, xml, force)
    fit_xml(xml, log, force)
    outdir = OUT / slug
    single_r = make_singlehit_figures(xml, outdir / f"kang24_{slug}_singlehit_delta", f"4TF full seed {seed}", force)
    fig4j_r = make_fig4j(xml, outdir / "fig4j", force)
    fig5_r = make_fig5(xml, outdir / "fig5", force)
    metrics = parse_log(log)
    return {
        "kind": "seed",
        "id": seed,
        "xml": str(xml),
        "log": str(log),
        "singlehit_r": single_r,
        "fig4j_multihit_r": fig4j_r,
        "fig5_at_r": fig5_r,
        **metrics,
    }


def run_scramble(rep: int, seed_for_annealer: int, force: bool = False) -> dict[str, object]:
    slug = f"4tf_scramble{rep}"
    template = WORK / f"kang24_{slug}_template.xml"
    scrambled = WORK / f"kang24_{slug}.xml"
    log = WORK / "logs" / f"{slug}.log"
    outdir = OUT / slug
    build_4tf_xml(seed_for_annealer, template, force)
    scramble_xml(template, scrambled, force)
    fit_xml(scrambled, log, force)
    # Evaluate fitted scrambled-label model against the true, unscrambled MPRA labels.
    single_r = make_singlehit_figures(scrambled, outdir / f"kang24_{slug}_true_labels_delta", f"4TF scramble {rep}", force)
    metrics = parse_log(log)
    return {
        "kind": "scramble",
        "id": rep,
        "annealer_seed": seed_for_annealer,
        "xml": str(scrambled),
        "log": str(log),
        "singlehit_r_against_true_labels": single_r,
        **metrics,
    }


def write_summary(rows: list[dict[str, object]]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    all_keys = sorted({key for row in rows for key in row})
    csv_path = OUT / "fig45_parallel_repeats_summary.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=all_keys)
        writer.writeheader()
        writer.writerows(rows)

    df = pd.DataFrame(rows)
    lines = ["# Fig. 4/5 4TF Repeats and Scramble Controls", "", "Competition mode: `false`; self-competition: `true`.", ""]
    for kind in ["seed", "scramble"]:
        sub = df[df["kind"].eq(kind)] if "kind" in df else pd.DataFrame()
        lines.append(f"## {kind}")
        lines.append("")
        lines.append(f"- Completed: `{len(sub)}`")
        for col in ["singlehit_r", "fig4j_multihit_r", "fig5_at_r", "singlehit_r_against_true_labels", "final_score", "runtime_seconds"]:
            if col in sub and sub[col].notna().any():
                lines.append(f"- {col}: mean `{sub[col].mean():.6g}`, median `{sub[col].median():.6g}`, min `{sub[col].min():.6g}`, max `{sub[col].max():.6g}`")
        lines.append("")
    lines.append(f"CSV: `{csv_path.name}`")
    (OUT / "fig45_parallel_repeats_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", type=int, default=5, help="Run 4TF fitting seeds 0..N-1.")
    parser.add_argument("--scrambles", type=int, default=10, help="Run scrambled-label controls 0..N-1.")
    parser.add_argument("--workers", type=int, default=min(6, os.cpu_count() or 4))
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--only", choices=["seeds", "scrambles"])
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    tasks = []
    if args.only in (None, "seeds"):
        for seed in range(args.seeds):
            tasks.append(("seed", seed))
    if args.only in (None, "scrambles"):
        for rep in range(args.scrambles):
            tasks.append(("scramble", rep))

    rows: list[dict[str, object]] = []
    print(f"tasks={len(tasks)} workers={args.workers}")
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        future_map = {}
        for kind, idx in tasks:
            if kind == "seed":
                future = executor.submit(run_seed, idx, args.force)
            else:
                future = executor.submit(run_scramble, idx, 1000 + idx, args.force)
            future_map[future] = (kind, idx)

        for future in concurrent.futures.as_completed(future_map):
            kind, idx = future_map[future]
            try:
                row = future.result()
                rows.append(row)
                write_summary(rows)
                print(f"DONE {kind} {idx}")
            except Exception as exc:
                rows.append({"kind": kind, "id": idx, "error": repr(exc)})
                write_summary(rows)
                print(f"FAILED {kind} {idx}: {exc}", file=sys.stderr)

    write_summary(rows)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
