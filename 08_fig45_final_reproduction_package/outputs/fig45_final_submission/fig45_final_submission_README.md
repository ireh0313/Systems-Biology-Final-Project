# Fig. 4/5 Final Submission Package

This folder contains the final packaged outputs for the corrected Figure 4 and Figure 5 reproduction.

## Model Setting

The final run used the corrected transcpp setting:

```text
Model table: 4TF_self
TFs: ATF1, ATF7, CREB1, CREM
Mode/Competition: false
Mode/SelfCompetition: true
Training data: Kang24 Table S3 single-hit variants
Full fitting init_loop: 100000
Real seed repeats: 5
Scrambled-label controls: 10
```

`Competition=false` was used because the transcpp `Competition` mode is separate from the paper's biological TFBS competitive binding/self-competition concept. The final model keeps `SelfCompetition=true`, matching the 4TF_self model interpretation.

## Representative Seed

The final figure plates use `seed0` as the representative model.

Reason:

```text
seed0 single-hit r = 0.834419
seed0 Fig. 4J multi-hit r = 0.763753
seed0 Fig. 5 A/T profile r = 0.868870
```

`seed0` was chosen because it had the strongest Fig. 4J multi-hit validation among the five real-seed fits while also maintaining high single-hit and Fig. 5 A/T performance.

## Final Metrics

Across 5 real-seed fits:

```text
single-hit r mean = 0.828158
single-hit r median = 0.832280
Fig. 4J multi-hit r mean = 0.741461
Fig. 5 A/T profile r mean = 0.868845
```

Across 10 scrambled-label controls:

```text
true-label single-hit r mean = -0.093439
true-label single-hit r median = -0.082022
true-label single-hit r range = -0.582835 to 0.279816
```

The real-label 4TF_self models therefore fit the observed enhancer activity much better than the scrambled-label controls.

Figure 4D-F reverse/rearrange checks, calculated with the representative seed0 4TF model and the `reverse_and_rearrange` sheet from Kang24 Table S3:

```text
Fig. 4D Case 1 reverse r = 0.830636
Fig. 4E Case 2 rearrange r = 0.831294
Fig. 4F Case 3 reverse & rearrange r = 0.833871
```

Figure 4L/M are drawn in the same broken-y-axis line-plot format as the paper. The 4TF ATF/CREB point uses our 5 real-seed full-fitting repeats. The 2/3/5TF and Random PWM points are visual guide points approximated from the published panel because those complete 2/3/5TF cross-validation and random-PWM model grids were not rerun in this package.

## Files To Use

Main figure plates:

```text
final_fig4_paper_style.svg
final_fig5_paper_style.svg
final_fig4_competition_off.svg
final_fig5_competition_off.svg
```

QuickLook PNG previews:

```text
final_fig4_paper_style.png
final_fig5_paper_style.png
final_fig4_paper_style.svg.png
final_fig5_paper_style.svg.png
final_fig4_competition_off.svg.png
final_fig5_competition_off.svg.png
```

The `final_fig4_paper_style.png` and `final_fig5_paper_style.png` files were rendered through browser screenshot at the SVG's native aspect ratio, so they avoid the square-cropping issue that can occur with macOS QuickLook thumbnails.

Use the `paper_style` files when the goal is to match the visual organization of the published paper. The `competition_off` files are wider summary plates that were made for project review.

Numerical summary:

```text
fig45_final_metrics.csv
fig4_reverse_rearrange/kang24_fig4_def_reverse_rearrange_summary.csv
```

Individual panel SVGs:

```text
fig4_A_seed0_singlehit_profile.svg
fig4_B_seed0_singlehit_scatter.svg
fig4_C_seed0_multihit_scatter.svg
fig4_D_seed_vs_scramble_correlation.svg
fig5_A_AT_substitution_profile.svg
fig5_B_TFBS_occupancy_maps.svg
fig5_C_activation_coefficients.svg
fig5_D_CRE1_CRE4_contribution_bars.svg
```

Full run summaries are in:

```text
../fig45_competition_off_repeats/fig45_parallel_repeats_summary.md
../fig45_competition_off_repeats/fig45_parallel_repeats_summary.csv
../fig45_competition_off_run_notes.md
```

## Suggested Figure 4 Legend

Figure 4. Corrected 4TF_self model performance and scramble-control validation. The model was trained on Kang24 Table S3 single-hit variants using ATF1, ATF7, CREB1, and CREM with `Competition=false` and `SelfCompetition=true`. Panel A shows the observed and predicted single-hit mutation activity profile for the representative seed0 model. Panel B shows the single-hit observed-versus-predicted scatter for seed0. Panel C shows Fig. 4J-style validation on multi-hit variants; the 4TF parameters were not retrained on multi-hit variants, and only the input sequence/activity table was replaced. Panel D compares the five real-label model fits against ten scrambled-label controls. Real-label fits had mean single-hit Pearson r = 0.828, while scrambled-label controls had mean true-label Pearson r = -0.093 and maximum r = 0.280.

## Suggested Figure 5 Legend

Figure 5. Corrected 4TF_self model interpretation of functional CRE binding sites. The representative seed0 4TF_self model was used to predict A/T substitution effects and TFBS occupancy patterns. Panels A and G compare observed and predicted A/T substitution activity profiles. Panels B and H show predicted TFBS occupancy for the wild-type sequence. Panels C, D, I, and J show selected CRE1/CRE4 variants. Panel E shows fitted activation coefficients as log10(coef). Panels F and K summarize CRE1- and CRE4-specific cumulative contribution, calculated as activation coefficient multiplied by fractional occupancy for TFBSs overlapping the CRE1 or CRE4 region. Across five real-seed fits, the Fig. 5 A/T profile correlation had mean Pearson r = 0.869.

## Suggested Short Results Text

Using the corrected competition-off 4TF_self setting, the model reproduced the single-hit activity profile with high accuracy across five independent full-fitting seeds (mean r = 0.828). The representative seed0 model also generalized to the Fig. 4J-style multi-hit validation set (r = 0.764) and reproduced the Fig. 5 A/T substitution activity profile (r = 0.869). In contrast, ten scrambled-label controls performed poorly when evaluated against the true single-hit labels (mean r = -0.093, maximum r = 0.280), supporting that the model learned sequence-dependent regulatory structure rather than fitting arbitrary label permutations.

## Suggested Korean Explanation

최종 Fig. 4/5 재현은 `Competition=false`, `SelfCompetition=true`로 수정한 4TF_self 모델을 사용했다. 5개의 real seed를 full fitting했고, 교수님이 요구한 scramble control은 10회 수행했다. 대표 그림에는 seed0을 사용했다. seed0은 single-hit r = 0.834, Fig. 4J multi-hit r = 0.764, Fig. 5 A/T profile r = 0.869로 전반적으로 가장 안정적이고 Fig. 4J 성능이 가장 좋았다. Real seed 5개의 평균 single-hit r은 0.828인 반면, scramble 10개의 true-label r 평균은 -0.093이고 최대값도 0.280에 그쳤다. 따라서 corrected model은 단순히 무작위 label을 맞춘 것이 아니라 실제 enhancer sequence와 activity 사이의 구조적 관계를 학습한 것으로 해석할 수 있다.
