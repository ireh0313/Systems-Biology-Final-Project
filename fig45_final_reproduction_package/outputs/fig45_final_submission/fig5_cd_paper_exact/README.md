# Figure 5C,D focused reproduction

This folder contains the focused reproduction of Kang et al. Figure 5C and Figure 5D.

## What is shown

- Figure 5C: CRE1 single A-substitution variants 11A, 14A, and 17A.
- Figure 5D: CRE4 single A-substitution variants 69A, 72A, and 75A.
- Each mini-panel shows predicted TF binding sites for ATF7, ATF1, CREM, and CREB1 using our corrected 4TF self-competition model.
- The bottom sequence track highlights the CRE-like motif in red and the substituted A position in yellow.

## Data source

- Model XML: `work/kang24/fig45_competition_off_repeats/kang24_4tf_seed0.xml`
- Substitution activity table: `outputs/fig45_competition_off_repeats/4tf_seed0/fig5/kang24_fig5_AT_substitution_table.csv`
- Binding-site calls exported here: `fig5_cd_variant_tfbs_occupancy_sites.csv`
- Summary metrics exported here: `fig5_cd_paper_exact_metrics.csv`

## Files

- `fig5_cd_paper_exact.svg`: editable vector version.
- `fig5_cd_paper_exact.png`: standard-resolution preview.
- `fig5_cd_paper_exact_4x.svg`: 4x-scale vector canvas.
- `fig5_cd_paper_exact_4x.png`: high-resolution PNG for reports/slides.
- `make_fig5_cd_paper_exact.py`: script used to generate the figure and CSV files.

