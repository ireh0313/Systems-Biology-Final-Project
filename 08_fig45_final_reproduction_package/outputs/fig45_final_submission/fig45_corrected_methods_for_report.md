# Corrected Methods For Reproducing Figures 4 and 5

We reproduced the Figure 4 and Figure 5 analyses using the Kang24 synthetic CRE enhancer model implemented in `transcpp`. The final analysis used the corrected 4TF_self model configuration with ATF1, ATF7, CREB1, and CREM. The model was trained on the single-hit variant measurements from Kang24 Table S3. We used full fitting with `init_loop=100000` and repeated the real-label fitting with five independent seeds.

The final XML model setting was:

```text
Mode/Competition = false
Mode/SelfCompetition = true
Mode/ScoreFunction = sse
Mode/ScaleData = false
Model = 4TF_self
TFs = ATF1, ATF7, CREB1, CREM
```

We set `Mode/Competition=false` because the transcpp `Competition` option is not the same as the biological TFBS competitive binding/self-competition mechanism described in the paper. The biological self-competition model was represented using `Mode/SelfCompetition=true`.

For Figure 4-style model-performance analysis, we first evaluated the trained 4TF_self model on the single-hit training variants. We plotted observed versus predicted activity values and calculated Pearson correlation between observed and predicted single-hit activity changes. We then performed a Fig. 4J-style multi-hit validation by keeping the fitted 4TF model parameters fixed and replacing only the input sequence/activity table with the multi-hit variants from Kang24 Table S3. No retraining was performed on the multi-hit variants.

For the scramble-control analysis, we generated ten scrambled-label controls using the `transcpp` scramble utility. Each scrambled control used the same model structure but with permuted `RateData` labels. After fitting to the scrambled labels, each fitted scrambled model was evaluated against the true, unscrambled single-hit labels. This allowed us to test whether the model could recover real sequence-dependent activity patterns rather than merely fitting arbitrary label assignments.

For Figure 5-style mechanistic interpretation, we used the representative seed0 4TF_self model. We predicted the effects of A and T substitutions across the CRE sequence and compared the predicted substitution profile with the observed profile. We also used transcpp `unfold`-derived occupancy outputs to visualize predicted TFBS occupancy around CRE1 and CRE4. We summarized a contribution proxy as activation coefficient multiplied by fractional occupancy. CRE1 and CRE4 were mapped to full-construct coordinates 31-38 and 89-96, corresponding to single-hit scan positions 11-18 and 69-76.

The representative seed0 model was selected for the final figure plates because it had strong performance across all Figure 4/5 analyses and the highest Fig. 4J multi-hit validation correlation among the five real-seed fits:

```text
seed0 single-hit Pearson r = 0.834419
seed0 Fig. 4J multi-hit Pearson r = 0.763753
seed0 Fig. 5 A/T profile Pearson r = 0.868870
```

Across five real-label fitting seeds, the single-hit mean Pearson r was 0.828158, the Fig. 4J multi-hit mean Pearson r was 0.741461, and the Fig. 5 A/T substitution profile mean Pearson r was 0.868845. Across ten scrambled-label controls, the mean Pearson r against the true single-hit labels was -0.093439, with a maximum of 0.279816. Therefore, the corrected 4TF_self model reproduced the real activity data much better than the scrambled-label controls.

The final packaged figures and numerical summaries are in:

```text
outputs/fig45_final_submission/
```

The full repeat-level summary is in:

```text
outputs/fig45_competition_off_repeats/fig45_parallel_repeats_summary.csv
outputs/fig45_competition_off_repeats/fig45_parallel_repeats_summary.md
```
