# Figure 4J Paper-Style Reproduction

This folder contains the focused Figure 4J reproduction in the same scatter-plot style as the Kang24 paper.

Files:

```text
fig4j_our_model_paper_exact.svg
fig4j_our_model_paper_exact.png
fig4j_our_model_paper_exact_4x.svg
fig4j_our_model_paper_exact_4x.png
fig4j_our_model_paper_exact_metrics.csv
make_fig4j_paper_exact.py
```

Run from the project root:

```bash
/Users/jisoochoi/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 work/kang24/make_fig4j_paper_exact.py
```

The plotted values use our reproduced 4TF_self multi-hit prediction output:

```text
observed x-axis: MPRA multi-hit data
predicted y-axis: Multi-hit prediction
```

The visual style follows the paper's Figure 4J format: gray panel background, white grid, blue points, panel letter J, axis labels, and Pearson R annotation.
