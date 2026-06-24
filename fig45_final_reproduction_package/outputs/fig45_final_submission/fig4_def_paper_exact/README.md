# Figure 4D-F Paper-Style Reproduction

This folder contains the focused Figure 4D-F reproduction in the same column-style format as the Kang24 paper.

Files:

```text
fig4_def_our_model_paper_exact.svg
fig4_def_our_model_paper_exact.png
fig4_def_our_model_paper_exact_metrics.csv
make_fig4_def_paper_exact.py
```

Run from the project root:

```bash
/Users/jisoochoi/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 work/kang24/make_fig4_def_paper_exact.py
```

The plotted values use the reproduced 4TF_self model outputs:

```text
Panel D: reverse sequence
Panel E: rearranged sequence
Panel F: reverse and rearranged sequence
```

The visual style follows the paper's Figure 4D-F panel format: vertical D/E/F layout, `Our model` title, gray panel background, white grid, inset scatter, and Pearson R annotation.
