# Figure 4L/M Meaning-Focused Reproduction

This figure is the recommended replacement for Figure 4L/M using only the data actually computed in this project.

Why this format:

```text
The paper's Fig. 4L/M plots 2, 3, 4, and 5 TF cross-validation curves plus Random PWM controls.
This project computed repeated full-fitting results for 4TF_self only, plus scrambled-label controls.
Therefore, a 2/3/4/5TF line plot would require inventing uncomputed values.
```

What this figure shows:

```text
Panel L: single-hit training-set Pearson R for 4TF_self real seeds versus scrambled-label controls.
Panel M: multi-hit validation-set Pearson R for 4TF_self real seeds versus scrambled-label controls.
```

Interpretation:

```text
The 4TF_self model gives consistently high correlations on both the training set and the multi-hit validation set, while scrambled-label controls stay near zero or negative. This preserves the core meaning of Fig. 4L/M: the real TF model captures sequence-dependent regulatory structure and generalizes beyond the fitting data.
```

Files:

```text
fig4lm_meaningful_4tf_vs_scramble.svg
fig4lm_meaningful_4tf_vs_scramble.png
fig4lm_meaningful_4tf_vs_scramble_4x.svg
fig4lm_meaningful_4tf_vs_scramble_4x.png
fig4lm_meaningful_4tf_vs_scramble_metrics.csv
make_fig4lm_meaningful.py
```
