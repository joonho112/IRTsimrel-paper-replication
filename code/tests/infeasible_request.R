# Reproduce one retained request whose target exceeds the detected range.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

source("code/00_setup.R")
source("code/01_simulation/validation.R")
use_irtsimrel("0.3.0"); load_scientific_helpers()
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
saved <- read_result("results", "s2", "calibration.csv")
reference <- saved[saved$eqc_status == "infeasible" & saved$pool == "parametric", ][1, ]
stopifnot(nrow(reference) == 1L, !is.na(reference$cond_id))
design <- read_result("results", "s2", "design_final.csv")
row <- design[design$cond_id == reference$cond_id, ]
context <- list(cfg = yaml::read_yaml(repo_path("config", "validation_design.yml")),
                pool_data = NULL)
value <- calibration_worker(row, context)
stopifnot(value$eqc_status == "infeasible",
          value$eqc_achieved < value$target_rho,
          abs(value$eqc_c_star - reference$eqc_c_star) < 1e-8,
          abs(value$eqc_achieved - reference$eqc_achieved) < 1e-8)
write_result(value, repo_path("output", "verification", "infeasible-request.csv"))
cat("Infeasible status, returned multiplier, and achieved index reproduced.\n")
