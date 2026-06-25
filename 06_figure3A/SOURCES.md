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

Unzip the archive into `data/` to restore the `data/rate_invert/` folder.

## Motif Sources

- CREB1 PWM: based on the original paper's supplementary PWM information.
- ATF/CREB family TF PWMs: based on the paper's supplementary motif information.
- Non-family TF PWMs: obtained from JASPAR.
- Homo sapiens motifs were prioritized because the original MPRA experiment used HEK293T cells.
- When exact human motifs were unavailable, human complex motifs or Mus musculus motifs were used as fallback choices.

## Reproduction Differences

The original paper trained more replicates per model group. This reproduction uses the available fitted outputs in this project folder, so the numerical values may not exactly match the published figure.
