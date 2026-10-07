# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

# Recompute published summaries and verify the distributed inputs.
# This command needs no calibration package and runs no response simulations.
source("code/00_setup.R")
require_packages("digest")
source(repo_path("code", "lib", "fct_reliability.R"))
checks <- list()
check <- function(name, value, expected, tolerance = 0) {
  passed <- length(value) == 1L && is.finite(value) &&
    abs(value - expected) <= tolerance
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name,
    observed = value, expected = expected, tolerance = tolerance, passed = passed
  )
  invisible(passed)
}
# Checksums detect an incomplete download or a changed retained input.
manifest_path <- repo_path("manifest", "files.csv")
if (!file.exists(manifest_path)) stop("Missing release manifest.")
m <- read.csv(manifest_path)
actual <- vapply(m$path, function(p) {
  full <- repo_path(p)
  if (!file.exists(full)) {
    return(NA_character_)
  }
  digest::digest(file = full, algo = "sha256")
}, character(1))
check("release files with matching SHA-256", sum(actual == m$sha256, na.rm = TRUE), nrow(m))

s2 <- read_result("results", "s2", "calibration.csv")
success <- s2[s2$eqc_status == "ok", ]
check("calibration requests", nrow(s2), 748)
check("unique calibration identifiers", length(unique(s2$cond_id)), 748)
check("successful calibrations", nrow(success), 743)
check("infeasible calibrations", sum(s2$eqc_status == "infeasible"), 5)
check(
  "EQC mean absolute solver residual", mean(abs(success$eqc_delta)),
  1.53762532914081e-9, 1e-20
)
check(
  "EQC maximum absolute solver residual", max(abs(success$eqc_delta)),
  3.80832176904988e-8, 1e-19
)
check("fixed-form package version", sum(s2$package_version == "0.3.0"), 748)
# The data table is independently recomputed from the individual node-study runs.
s3 <- read_result("results", "s3", "m_sensitivity.csv")
check("node-study runs", nrow(s3), 1680)
summary <- read_result("tables", "T_s3_m.csv")
for (i in seq_len(nrow(summary))) {
  z <- s3[s3$M == summary$M[i] & s3$status == "ok", ]
  check(
    paste("fresh-draw MAE at M", summary$M[i]), mean(abs(z$holdout_error)),
    summary$holdout_mae[i], 1e-14
  )
  check(paste("node-study runs at M", summary$M[i]), nrow(z), 240)
}
sp <- read_result("results", "s2", "superpopulation_v031.csv")
check("item-population requests", nrow(sp), 380)
check("item-population package version", sum(sp$sac_superpop_package_version == "0.3.1"), 380)
# Arm-A rates use dataset-level model rows, not the displayed rounded percentages.
a <- read_result("results", "s6", "arm_a_runs.csv")
at <- read_result("tables", "T_s6_arm_a.csv")
keys <- c("target_rho", "te_mean", "te_sd", "te_scale", "method")
for (i in seq_len(nrow(at))) {
  keep <- rep(TRUE, nrow(a))
  for (key in keys) keep <- keep & a[[key]] == at[[key]][i]
  z <- a[keep & !is.na(a$estimate), ]
  check(paste("Arm-A rejection rate", i), mean(z$p.value < .05), at$reject_rate[i], 1e-12)
}
# Mathematical properties catch changes that aggregate-file comparisons cannot.
source(repo_path("code", "tests", "numerical_checks.R"))
result <- do.call(rbind, checks)
write_result(result, repo_path("output", "verification", "checks.csv"))
cat(sprintf("Verification: %d / %d checks passed.\n", sum(result$passed), nrow(result)))
if (!all(result$passed)) {
  print(result[!result$passed, ], row.names = FALSE)
  stop("Verification failed. Inspect output/verification/checks.csv.")
}
