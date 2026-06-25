# Figure 3A Reproduction

This folder contains code and intermediate `unfold --rate --invert` outputs used to reproduce Figure 3A from Kang and Kim, iScience 2024.

## Purpose

Figure 3A compares a CREB1 self-competition baseline model with 2TF models in which one additional TF is added to CREB1. The goal is to test whether non-family TFs or ATF/CREB family TFs improve the model's ability to explain single-hit mutation activities of the CRE enhancer.

## Folder Structure

```text
06_figure3A/
├─ data/
│  ├─ Figure3A_excluded_nan_runs.txt
│  ├─ mmc_files/
├─ outputs/
├─ scripts/
│  ├─ figure_metrics.py
│  └─ plot_figure3a.py
├─ data_rate_invert.zip
├─ README.md
├─ requirements.txt
└─ SOURCES.md
```

## Required External Input

Place the supplementary Table S3 Excel file in:

```text
data/mmc_files/1-s2.0-S2589004223028249-mmc4.xlsx
```

This file is not included here. Download it from the supplementary information of Kang and Kim, iScience 2024.

The intermediate `unfold --rate --invert` outputs are compressed in:

```text
data_rate_invert.zip
```

Before rerunning the plotting script, unzip this file into the `data/` folder so that `data/rate_invert/` is restored.

## How to Run

From this folder:

```bash
unzip data_rate_invert.zip -d data/
python scripts/plot_figure3a.py
```

Optional:

```bash
python scripts/plot_figure3a.py --show
python scripts/plot_figure3a.py --ymin 0.4 --ymax 0.9
```

## Outputs

The script writes:

```text
outputs/Figure3A_TF_addition.png
outputs/Figure3A_TF_addition.pdf
outputs/Figure3A_TF_addition_metrics.csv
```

## Notes

- CREB1-only self-competition models are used as the baseline.
- Family TFs use ATF/CREB family motifs.
- Non-family TFs were generated using JASPAR-derived motifs.
- Runs that ended with `s = -nan` are plotted as red X marks but are excluded from the boxplot statistics.
- Pearson's R is calculated between WT-normalized log2 experimental activity and WT-normalized log2 model-predicted rate.
