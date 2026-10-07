<!-- Author: JoonHo Lee (jlee296@ua.edu) -->

# Reading the result files

| Field or group | Meaning |
|---|---|
| `cond_id` | Stable condition identifier; use it to match a rerun to the retained study |
| `seed`, `seed_id` | Actual design seed, or index used in the documented seed formula |
| `model` | `rasch`, `2pl`, `3pl_g20` (guessing .20), or `3pl_beta` (Beta(5,17) guessing) |
| `latent_shape`, `shape` | Normal, bimodal, positively skewed, or standardized heavy-tailed t5 distribution |
| `pool` | Parametric difficulties or empirical IRW difficulties |
| `n_items`, `n_persons`, `n_per_group` | Test length and person sample size; group sizes are per group |
| `target_rho` | Requested mean-information index; not an estimated score reliability |
| `eqc_c_star`, `c_star` | Common multiplier applied to baseline item discriminations |
| `eqc_achieved`, `rho_tilde` | Mean-information index evaluated using the specified nodes and variance |
| `eqc_delta`, `calibration_residual` | Achieved index minus target on the solver's own nodes |
| `holdout_error` | Index minus target on independent evaluation draws |
| `rho_psd` | Historical column name for the pointwise-information index |
| `w_bar` | Inverse-information index; missing when the population reciprocal expectation diverges |
| `w_bar_truncated` | Finite-node reciprocal-information calculation; not a population value in the divergent-tail case |
| `theta_var` | Variance used to define the index; inspect its source as well as its value |
| `sac_superpop_rho_at_expected_information` | Primary item-population index: transform of the mean information across forms |
| `sac_superpop_mean_random_form_reliability` | Descriptive mean of the separately transformed form indices |
| `sac_superpop_jensen_gap` | Difference between the preceding two quantities |
| `status`, `eqc_status` | Completion or feasibility status; do not remove failed rows before defining a denominator |
| `elapsed_sec`, `*_elapsed_sec` | Recorded runtime; expected to differ on another computer |
| `te_mean`, `te_sd` | Mean and across-item SD of the treatment-effect distribution under its stated scale convention |
| `te_scale` | `theta_fixed` holds latent-coordinate effects fixed; `logit_fixed` holds logit coefficients fixed |
| `method` | Constant or varying item-effect model fitted to the same response dataset |
| `estimate`, `std.error`, `p.value` | Fitted average treatment effect, its model SE, and the test p-value |
| `converged` | Historical fit flag; it is not the later Arm-B admissibility decision |
| `form_numerator`, `form_denominator` | Form-level rejection count and evaluable-item count for MH inference |
| `rep`, `draw` | Independent response replication or item-form draw within the stated design |

`data/precomputed/results/s2/design_final.csv` contains 748 accepted requests. Its predecessor, `feasibility.csv`, has 128 structural cells. The screen and calibration use separately seeded forms, so passing the screen does not guarantee every later form will reach the requested target.

`data/precomputed/results/s3/m_sensitivity.csv` has 1,680 rows: 12 structural cells, seven node counts, and 20 seeds. Its fresh-draw sample has 200,000 nodes in every condition.

`data/precomputed/results/s6/arm_a_runs.csv` has two model rows per generated response dataset. Counting both models as independent replications would give the wrong denominator.
