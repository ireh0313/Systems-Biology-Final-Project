# CRE enhancer transcpp reproduction notes

This repository contains the code and command notes used to reproduce selected modeling and visualization steps from the CRE enhancer study based on DOI `10.1016/j.isci.2023.108747`.

The code is organized by analysis step:

1. `01_build_neoparsa_transcpp`: build `neoParSA` and `transcpp` on Ubuntu/WSL.
2. `02_kim2013_example_check`: run the bundled Kim2013 example to verify the build.
3. `03_figure2A_2B`: reproduce the CRE sequence schematic and single-hit MPRA mutation-effect plots.
4. `04_figure2C`: reproduce CREB1-only baseline and 7-family model prediction plots.
5. `05_figure2D`: reproduce 7-family PWM sequence logos.
6. `common`: shared data-preparation scripts.

## Data Not Included

Large result files, fitted XML outputs, and the paper supplementary Excel files are not committed here by default. Place the supplementary files under:

```text
~/cre_repro/data/
```

Expected files:

```text
1-s2.0-S2589004223028249-mmc2.xlsx
1-s2.0-S2589004223028249-mmc3.xlsx
1-s2.0-S2589004223028249-mmc4.xlsx
```

The scripts assume the working directory layout used during the project:

```text
~/neoParSA-master
~/transcpp-master
~/cre_repro
~/cre_repro/data
~/cre_repro/work
~/cre_repro/results
~/cre_repro/scripts
```

## Python Environment

Create a virtual environment before running the plotting/data scripts:

```bash
cd ~/cre_repro
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -r /path/to/this/repo/requirements.txt
```

## Basic Reproduction Order

```bash
# 1. Build the old C++ tools.
bash 01_build_neoparsa_transcpp/build_neoparsa_transcpp_ubuntu.sh

# 2. Verify the executable with Kim2013.
bash 02_kim2013_example_check/run_kim2013_smoke_test.sh

# 3. Prepare supplementary data tables.
python common/export_tables.py
python common/prepare_single_hit.py

# 4. Run Figure 2A/2B scripts.
python 03_figure2A_2B/scan_cre_sites_pwm.py
python 03_figure2A_2B/reproduce_fig2_paper_method.py
python 03_figure2A_2B/plot_fig2a_complete.py
python 03_figure2A_2B/plot_fig2b_position_margin_grid_final.py

# 5. Run Figure 2C scripts after fitted model XML/rate files are available.
python 04_figure2C/plot_fig2c_baseline_clean_inset.py
python 04_figure2C/plot_fig2c_7tf_family_model_grid19_refined.py

# 6. Run Figure 2D sequence logos.
python 05_figure2D/plot_fig2d_sequence_logos.py
```

## Notes On Interpretation

Figure 2B uses measured MPRA expression values from Table S3 (`mmc4`). The plotted value is:

```text
log2(variant expression / WT expression)
```

Figure 2C uses both Table S2 (`mmc3`, PWM/PSSM values) and Table S3 (`mmc4`, sequences and measured expression). PWM values define TF binding scores in the model; measured expression values are the fitting target and evaluation reference.

Figure 2D uses only Table S2 (`mmc3`) to draw PWM sequence logos. It does not use MPRA expression values and does not require model fitting.

