# Figure 3B and Figure 3D Reproduction

This folder contains code and intermediate `unfold --rate --invert` outputs used to reproduce Figures 3B and 3D from Kang and Kim, iScience 2024.

## Purpose

Figures 3B and 3D evaluate how the number of ATF/CREB family TFs affects model performance in single-hit mutation fitting.

- Figure 3B: single-hit fitting without self-competition
- Figure 3D: single-hit fitting with self-competition

Both figures use the same single-hit MPRA data and the same seven ATF/CREB family TFs. The key difference is whether self-competition is included.

## Folder Structure

```text
07_figure3B_3D/
├─ data/
│  └─ mmc_files/
├─ outputs/
├─ scripts/
│  ├─ figure_metrics.py
│  └─ plot_figure3b_3d.py
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

Before rerunning the plotting script, unzip this file into the `data/` folder so that `data/figure3b_selfcompetition_off/` and `data/figure3d_selfcompetition_on/` are restored.

## How to Run

Generate both figures:

```bash
unzip data_rate_invert.zip -d data/
python scripts/plot_figure3b_3d.py
```

Generate one figure:

```bash
python scripts/plot_figure3b_3d.py --figure 3B
python scripts/plot_figure3b_3d.py --figure 3D
```

Show the plot window:

```bash
python scripts/plot_figure3b_3d.py --show
```

## Outputs

```text
outputs/Figure3B_selfcompetition_off.png
outputs/Figure3B_selfcompetition_off.pdf
outputs/Figure3B_selfcompetition_off_metrics.csv

outputs/Figure3D_selfcompetition_on.png
outputs/Figure3D_selfcompetition_on.pdf
outputs/Figure3D_selfcompetition_on_metrics.csv
```

## Notes

- The seven TFs are CREB1, CREB3, CREB5, CREM, ATF1, ATF4, and ATF7.
- All possible TF combinations are evaluated: 7, 21, 35, 35, 21, 7, and 1 combinations for 1TF through 7TF.
- The reproduction uses one fitted result per TF combination.
- Pearson's R is calculated between WT-normalized log2 experimental activity and WT-normalized log2 model-predicted rate.
