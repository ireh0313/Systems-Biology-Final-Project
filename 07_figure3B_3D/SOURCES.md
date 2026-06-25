# Sources and Data Provenance

## Original Paper

Kang, C.-K. and Kim, A.-R.  
Deep molecular learning of transcriptional control of a synthetic CRE enhancer and its variants.  
iScience 27, 108747, 2024.

## MPRA Activity Data

Single-hit mutation expression values are from supplementary Table S3:

```text
1-s2.0-S2589004223028249-mmc4.xlsx
```

The script uses the `single-hit` sheet and the `expression level` column.

## Model Output Files

The compressed archive `data_rate_invert.zip` contains files generated from fitted transcpp XML models using:

```bash
./unfold -i model.xml -s Output --rate --invert > model_rate_invert.txt
```

Unzip the archive into `data/` to restore `data/figure3b_selfcompetition_off/` and `data/figure3d_selfcompetition_on/`.

## Figure 3B

Figure 3B uses self-competition-off single-hit fitting outputs.

## Figure 3D

Figure 3D uses self-competition-on single-hit fitting outputs from `Fig3D_unfold_rates`.

The `2TF_full_partial` files were not used for the final Figure 3D reproduction because they are an incomplete partial set. The plotted Figure 3D uses:

- `1TF_full`
- `2TF_medium`
- `3TF_medium`
- `4TF_medium`
- `5TF_medium`
- `6TF_full`
- `7TF_full`

## Reproduction Differences

The original paper trained at least eight models for each TF combination. This reproduction uses one fitted result per TF combination, so exact numerical distributions may differ from the published figure.
