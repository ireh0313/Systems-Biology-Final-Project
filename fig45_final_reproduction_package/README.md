# SysBio Final Project: Competition-Off Final Reproduction Package

This is the corrected GitHub upload package based on the final DOCX report and slide PDF.

The selection rule for this folder is:

- Include the final `competition=off` reproduction workflow and its supporting files.
- Include all real-seed and scramble-control outputs described in the final report process.
- Include the final `paper_exact_4x` reproduced figure panels used in the report/slides.
- Exclude older or inconsistent reproduction branches such as `fullfit_final`, `fig4j_multihit`, and `fig5_analysis`.

## Main Contents

- `docs/`
  - Final report DOCX
  - Final slide PDF
- `outputs/fig45_competition_off_repeats/`
  - Corrected final workflow outputs
  - Includes `4tf_seed0` through `4tf_seed4`
  - Includes `4tf_scramble0` through `4tf_scramble9`
  - Includes `fig45_parallel_repeats_summary.csv/.md`
- `outputs/fig45_final_submission/`
  - Final figure panels and final submission notes
  - Includes the `paper_exact` and `paper_exact_4x` figure folders:
    - `fig4_def_paper_exact`
    - `fig4j_paper_exact`
    - `fig4lm_meaningful`
    - `fig4lm_paper_exact`
    - `fig5_ag_paper_exact`
    - `fig5_bh_paper_exact`
    - `fig5_cd_paper_exact`
    - `fig5_ij_paper_exact`
    - `fig5e_paper_exact`
    - `fig5_fk_paper_exact`
- `work/kang24/`
  - Scripts used for the corrected workflow and final paper-exact figure generation
  - Supplementary input files
  - Corrected `competition=off` model XML files for seeds/scrambles
  - Selected paper crop references used by paper-exact figure scripts

## Notes

- The large external TransCPP repository/build folder is not included. The scripts expect it at `work/repos/transcpp` if the full workflow is rerun.
- Run logs and `.prolix` fitting traces were excluded to avoid uploading bulky execution noise; the fitted XMLs, CSV summaries, and final figure outputs are included.
- The older broad folders `outputs/fullfit_final`, `outputs/fig4j_multihit`, and `outputs/fig5_analysis` are intentionally excluded because they do not match the final `competition=off` reproduced figures.
