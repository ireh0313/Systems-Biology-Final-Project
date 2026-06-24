# 04. Figure 2C

Figure 2C compares model-predicted activity changes with measured Figure 2B MPRA activity changes.

## Baseline Model

Baseline model:

- CREB1-only model
- Self-competition OFF
- Full annealing
- Uses Table S2 CREB1 PWM/PSSM for binding scores
- Uses Table S3 single-hit expression data as fitting/evaluation data

Input used during this project:

```text
~/CRE_1TF_selfoff_singlehit_fullanneal_CREB1_scramble_rate_invert.txt
```

Plot:

```bash
python plot_fig2c_baseline_clean_inset.py
```

Final baseline panel is the small inset labeled `B.M`.

## 7TF Family Model

7TF model:

- CREB1, CREB3, CREB5, CREM, ATF1, ATF4, ATF7
- Self-competition ON
- Competition OFF
- Full fitting performed on AWS

Input XML used during this project:

```text
~/cre_repro/results/aws_7TF_full/CRE_7TF_full_seed870001.xml
```

The script extracts predicted rates using:

```bash
~/transcpp-master/unfold -i ~/cre_repro/results/aws_7TF_full/CRE_7TF_full_seed870001.xml --rate --invert
```

Predicted activity change:

```text
log2(mutant predicted rate / WT predicted rate)
```

Observed activity change:

```text
log2(mutant measured expression / WT measured expression)
```

Plot:

```bash
python plot_fig2c_7tf_family_model_grid19_refined.py
```

The reported `R`, `R^2`, and `RMSE` values compare the 7TF predicted log2 activity change against the measured Figure 2B log2 activity change for each substituted base.

