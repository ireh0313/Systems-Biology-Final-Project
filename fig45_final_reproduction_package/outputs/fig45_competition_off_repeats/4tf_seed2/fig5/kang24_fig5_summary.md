# Kang24 Figure 5-Style Functional Binding Site Analysis

This analysis uses the completed 4TF full fitted model and `unfold` outputs for sites and fractional occupancy.

- 4TF model A/T substitution Pearson r: `0.865157`
- Contribution proxy used here: `activation coefficient x fractional occupancy`.
- CRE1 and CRE4 coordinates are mapped from the single-hit scan: CRE1 `31-38`, CRE4 `89-96` in the full construct, corresponding to scan positions `11-18` and `69-76`.

## Activation coefficients

- ATF1: `1.04237`
- ATF7: `3.44698e-06`
- CREB1: `4.68702`
- CREM: `19.9195`

## Key variants

| Variant | Observed delta | Predicted delta |
|---|---:|---:|
| synCRE_Promega_0 | `0.000000` | `0.000000` |
| scanmut_single_pos_14_A | `-0.304006` | `-0.484986` |
| scanmut_single_pos_14_T | `-0.377070` | `-0.408170` |
| scanmut_single_pos_72_A | `-1.888969` | `-1.324058` |
| scanmut_single_pos_72_T | `-1.434403` | `-1.041166` |

## Outputs

- `kang24_fig5_AT_substitution_profile.svg`
- `kang24_fig5_TFBS_occupancy_maps.svg`
- `kang24_fig5_activation_coefficients.svg`
- `kang24_fig5_CRE1_CRE4_contribution_bars.svg`
- `kang24_fig5_sites_occupancy_contribution.csv`
- `kang24_fig5_AT_substitution_table.csv`
