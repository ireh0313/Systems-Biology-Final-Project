# 03. Figure 2A and 2B

Figure 2A shows the synthetic CRE enhancer sequence, known CRE1-4 sites, and the cryptic region.

Figure 2B shows measured single-hit mutation effects by substituted base.

## Input Data

The main input is Table S3 from the paper supplement:

```text
~/cre_repro/data/1-s2.0-S2589004223028249-mmc4.xlsx
```

The PWM/PSSM supplement is used only as motif-scan support:

```text
~/cre_repro/data/1-s2.0-S2589004223028249-mmc3.xlsx
```

## Data Preparation

```bash
cd ~/cre_repro
source .venv/bin/activate

python /path/to/repo/common/export_tables.py
python /path/to/repo/common/prepare_single_hit.py
```

Important processing choices:

- WT row: `synCRE_Promega_0`
- Full WT sequence length: 197 bp
- CRE scan region: 87 bp, extracted from 0-based full-sequence slice `20:107`
- Public Table S3 contains WT plus 260 single-hit variants.
- The expected `scanmut_single_pos_1_T` row is absent from the public supplement and was not synthesized.

## Final Figure Commands

```bash
python scan_cre_sites_pwm.py
python reproduce_fig2_paper_method.py
python plot_fig2a_complete.py
python plot_fig2b_position_margin_grid_final.py
```

## Final Cryptic Region Rule

The final paper-like cryptic region was selected from measured MPRA effects:

```text
outside canonical CRE sites
all three substitutions have positive log2 activity change
mean log2 activity change >= 0.1
standard deviation across substitutions <= 0.15
longest contiguous region
```

Final region:

```text
scan positions 63-66, sequence CCCC
```

The plotted mutation-effect value is:

```text
log2(variant expression / WT expression)
```

not `log2(variant expression - WT expression)`.

