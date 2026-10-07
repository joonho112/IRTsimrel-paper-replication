<!-- Author: JoonHo Lee (jlee296@ua.edu) -->

# Release checks

Checked on 7 October 2026 with R 4.6.0 on macOS. The dependency record is in `environment/packages.csv` and `renv.lock`.

| Check | Executed scope | Result |
|---|---|---|
| Retained evidence and numerical identities | 179 checks, including file hashes, EQC residuals, node-study summaries, Arm-A rates, invariance, tail handling, and invalid inputs | Passed |
| Feasibility screen | All 128 structural cells, both difficulty pools | Recovered all 748 accepted requests with identical identifiers and seeds; maximum range discrepancy below 4.0e-15 |
| Fixed-form calibration | 32 smoke conditions, both pools | All 96 compared values within 1e-8 of retained values |
| Infeasible request | One additional parametric request at target .95 | Infeasible status, best detected multiplier, and achieved index reproduced |
| Node-count study | Four smoke conditions | All 12 compared values within 1e-8 |
| Fixed-form algorithm comparison | Two smoke conditions | All six compared values within 1e-8 |
| Stochastic-step study | Two smoke conditions | All eight compared values within 1e-8 |
| Item-population calibration | One smoke condition, version 0.3.1 | Both compared values within 1e-8 |
| Repeated-person sampling | One condition, five person sample sizes, ten repetitions each | Completed; calibrated bank identity checked |
| Parameter recovery | Two conditions, two response replications each | Completed with TAM |
| Item recipe atlas | Four shape conditions, two forms each | Completed |
| Treatment-effect Arm A | Four datasets, eight fitted-model rows | Completed |
| Treatment-effect Arm B | Two datasets, four initial fitted-model rows | Completed; later retry analysis not rerun |
| MH form-cluster study | One structural cell, two forms, 96 retained aggregate rows | Completed |
| Score comparison | Four conditions, two response replications each | Completed |
| Checkpoints and package isolation | Reuse, changed-context rejection, unsafe identifier rejection, wrong loaded-version rejection | Passed |
| Numerical figures | 18 PDF/PNG figures rebuilt | All 18 PNGs pixel-identical to the v6.3 manuscript figures on this machine |
| Numerical tables | 12 LaTeX files recalculated | Completed; final reference wrappers retained separately |
| Reader guide | Seven HTML chapters | Rendered successfully; all local image and asset links resolve |
| Relocated source ZIP | Extracted outside the original workspace; verification, frozen-package installation, worked example, 16-condition calibration, checkpoint checks, and full quick route | Passed; all 48 calibration comparisons within 1e-8 |
| Public IRW retrieval | Downloaded the pinned public source archive and extracted its difficulty data | Matched the frozen SHA-256; 6,143 rows and 145 instruments |
| Source preparation | Code parsing, authored-file headers, release paths, archive contents, and file-level authorship | Passed |

The eight scientific helper modules preserve their original computations. Their parsed expressions agree with the retained implementation after accounting for single-statement braces introduced by formatting. Package archives and retained numerical inputs preserve their original bytes.

The guide's HTML was checked through rendering and static asset validation; a browser appearance review was not performed. GitHub Actions is configured but has not been run on a remote repository. The complete historical Monte Carlo workload, a fresh cross-platform dependency restore, and the retained-only supplementary analyses were not rerun for this release. See `docs/reproduction-scope.md` for the distinction between executable studies, recalculated displays, and retained evidence.

To repeat the fast checks:

```sh
Rscript run_all.R verify
Rscript code/tests/checkpoint_checks.R
Rscript code/tests/infeasible_request.R
```

## Public release review

The public snapshot retains the study inputs, source packages, scientific outputs used by the manuscript, and reader-facing verification documentation. It omits runtime libraries, generated execution outputs, local working material, and 114 unused duplicate or working-record files from the preparation copy.

The public IRW difficulty export is included with its unchanged source license and a verified checksum. A fresh 16-condition IRW calibration run reproduced all 48 compared values within 1e-8. The reduced display-input set rebuilt all 18 numerical figures and 12 computational tables. The two supplementary person-sampling and treatment-effect figure entries in the reproduction map were checked against the manuscript and corrected.

The public snapshot passed all 179 verification checks. Its Git history starts with the reviewed release rather than retaining preparation records. The replication-package author is JoonHo Lee (jlee296@ua.edu); the original IRW attribution is preserved separately.
