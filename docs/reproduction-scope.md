<!-- Author: JoonHo Lee (jlee296@ua.edu) -->

# Reproduction scope

We distinguish three operations throughout the repository.

1. **Recalculate a display:** read the retained results and compute a figure, table, or verification statistic again. This does not regenerate responses.
2. **Rerun a design:** generate new numerical output with the paper's seeds and settings. The validation commands implement this operation in smoke and full modes.
3. **Inspect retained evidence:** read a result whose complete historical pipeline is not part of this release.

## Validation

The feasibility, fixed-form calibration, solver comparison, item-population calibration, integration-node, stochastic-step, repeated-person-sampling, and parameter-recovery designs have executable full paths. Feasibility reconstructs the accepted target requests and checks their identifiers and seeds against the retained design. Downstream commands deliberately start from that retained, versioned design so they can also be run independently.

The fixed-form results used IRTsimrel 0.3.0. The worked example and item-population study used 0.3.1. The repository installs these versions in separate project libraries and checks the loaded namespace. It does not replace either with a currently installed global version.

The paper's settings and seeds are preserved. Small floating-point differences can occur across R, compiler, BLAS, and model-fitting environments. Calibration comparisons use an absolute tolerance of 1e-8. This is a replay tolerance, not a new criterion for scientific success. Solver accuracy and independent-draw discrepancy are reported separately.

## Supplementary studies

- `recipes` regenerates the item-form atlas draws. The broader descriptive comparison with the literature's reported reliability corpus is retained-only; that corpus is not redistributed.
- `treatment` regenerates Arm A, including the primary latent-coordinate effect convention and the fixed-logit sensitivity. It uses the paper's independent target-specific seeds. The command does not first rerun the external Gilbert–Kim–Miratrix fidelity comparison; the paper's 72-comparison report is included for inspection.
- `treatment-b` regenerates the original Arm-B fits. It does **not** implement the later optimizer retry battery, fit-admissibility decisions, or boundary-aware publication pipeline. Those resolved model rows and sensitivity tables are retained inputs for the paper's supplementary displays. Initial fits must not be substituted for the resolved analysis.
- `dif` regenerates the 100-form-per-cell MH study, preserving paired streams and form-level rejection numerators and denominators. Historical 1,000/500-form point estimates and their logistic-regression summaries are retained-only. The bootstrap publication tables are also retained; the library exposes the form-cluster calculations for further work.
- `scores` reruns the score-comparison study with TAM. Smoke mode uses two response replications, full mode 200.
- The detailed theory grids, equalization-cost sweep, and external shape/reliability benchmarks are retained evidence. The verification command independently checks selected mathematical identities; it does not rerun the complete historical theory census.

## Displays

The quick route rebuilds every numerical figure used in the manuscript and supplement, and 12 selected numerical LaTeX tables from the distributed rows or summaries. Some plots read retained aggregate tables; they are display reproductions, not independent reruns of the simulation or of every intermediate aggregation.

The remaining table fragments and the conceptual forward/inverse diagram are supplied as reference source. Table captions, labels, and layout can differ between a computational table builder and the final typeset fragment. The reproduction map states which role each file plays.

The manuscript itself is not duplicated here. This repository accompanies v6.3; it is not a new manuscript revision.
