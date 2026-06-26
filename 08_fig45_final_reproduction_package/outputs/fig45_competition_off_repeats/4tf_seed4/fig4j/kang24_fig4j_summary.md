# Kang24 Figure 4J Multi-Hit Prediction

This reproduces the Figure 4J-style multi-hit validation using the completed 4TF full fitted model.

- Model XML with multi-hit genes: `kang24_4tf_full_seed0_multihit.xml`
- Multi-hit variants included: `722`
- Pearson r: `0.710158`
- RMSE: `0.935265`

The fitted 4TF parameters were not retrained on multi-hit data. Only the gene sequence table and observed RateData were replaced with the `multi-hit` sheet from Kang24 Table S3.

## Outputs

- `kang24_fig4j_4tf_multihit_scatter.svg`
- `kang24_fig4j_4tf_multihit_observed_predicted.csv`
- `kang24_4tf_full_seed0_multihit.xml`
