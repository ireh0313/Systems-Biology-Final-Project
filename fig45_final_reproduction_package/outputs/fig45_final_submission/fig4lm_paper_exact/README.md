# Figure 4L/M Paper-Style Reproduction

This folder contains focused Figure 4L/M panels in the same broken-axis style as the Kang24 paper.

Files:

```text
fig4lm_our_values_paper_exact.svg
fig4lm_our_values_paper_exact.png
fig4lm_our_values_paper_exact_4x.svg
fig4lm_our_values_paper_exact_4x.png
fig4lm_our_values_paper_exact_metrics.csv
make_fig4lm_paper_exact.py
```

Run from the project root:

```bash
/Users/jisoochoi/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 work/kang24/make_fig4lm_paper_exact.py
```

Data used:

```text
Panel L real series: 5 real-seed 4TF_self single-hit training correlations
Panel L control series: 10 scrambled-label single-hit true-label correlations
Panel M real series: 5 real-seed 4TF_self multi-hit validation correlations
Panel M control series: 10 scrambled-label models evaluated on true multi-hit labels
```

Important note:

```text
The cyan series is labeled scrambled-label controls, not Random PWMs, because Random PWM model grids were not recomputed in this project.
```
