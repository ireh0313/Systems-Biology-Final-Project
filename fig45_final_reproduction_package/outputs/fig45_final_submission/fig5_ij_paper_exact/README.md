# Figure 5I,J focused reproduction

This folder contains the focused reproduction of Kang et al. Figure 5I and Figure 5J.

## What is shown

- Figure 5I: CRE1 single T-substitution variants 12T, 15T, and 18T.
- Figure 5J: CRE4 single T-substitution variants 70T, 73T, and 76T.
- Each mini-panel shows predicted TF binding sites for ATF7, ATF1, CREM, and CREB1 using our corrected 4TF self-competition model.
- The bottom sequence track highlights the CRE-like motif in red, the substituted T position in yellow, and the CRE4/cryptic region marker in blue for panel J.

## Data source

- Model XML: `work/kang24/fig45_competition_off_repeats/kang24_4tf_seed0.xml`
- Substitution activity table: `outputs/fig45_competition_off_repeats/4tf_seed0/fig5/kang24_fig5_AT_substitution_table.csv`
- Binding-site calls exported here: `fig5_ij_variant_tfbs_occupancy_sites.csv`
- Summary metrics exported here: `fig5_ij_paper_exact_metrics.csv`

## Files

- `fig5_ij_paper_exact.svg`: editable vector version.
- `fig5_ij_paper_exact.png`: standard-resolution preview.
- `fig5_ij_paper_exact_4x.svg`: 4x-scale vector canvas.
- `fig5_ij_paper_exact_4x.png`: high-resolution PNG for reports/slides.
- `make_fig5_ij_paper_exact.py`: script used to generate the figure and CSV files.

