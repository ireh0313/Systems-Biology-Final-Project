# Figure 5F,K focused reproduction

This folder contains the focused reproduction of Kang et al. Figure 5F and Figure 5K.

## What is shown

- Figure 5F: cumulative CRE1/CRE4 contribution after A-substitution variants.
- Figure 5K: cumulative CRE1/CRE4 contribution after T-substitution variants.
- Bars are stacked by individual TFBS contribution, grouped in the order CREB1, CREM, ATF1, and ATF7.
- The contribution value is calculated as activation coefficient multiplied by fractional TFBS occupancy.
- Grey headers separate CRE1 and CRE4, and black connector lines link cumulative site-level contribution boundaries between variants.

## Data source

- Model XML: `work/kang24/fig45_competition_off_repeats/kang24_4tf_seed0.xml`
- Site-level contribution table exported here: `fig5_fk_site_level_contributions.csv`
- Stacked bar values exported here: `fig5_fk_stacked_bar_values.csv`

## Files

- `fig5_fk_paper_exact.svg`: editable vector version.
- `fig5_fk_paper_exact.png`: standard-resolution preview.
- `fig5_fk_paper_exact_4x.svg`: 4x-scale vector canvas.
- `fig5_fk_paper_exact_4x.png`: high-resolution PNG for reports/slides.
- `make_fig5_fk_paper_exact.py`: script used to generate the figure and CSV files.
