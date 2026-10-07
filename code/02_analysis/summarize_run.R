# Summarize a newly executed study without altering the paper's retained inputs.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

summarize_run <- function(study, mode) {
  path <- repo_path("output", "runs", study, mode, "results.csv")
  if (!file.exists(path)) stop("Run the study first: ", study, " (", mode, ").")
  x <- read.csv(path)
  if (study == "calibration") {
    out <- do.call(rbind, lapply(split(x, x$eqc_status), function(z) {
      data.frame(
        status = z$eqc_status[1], requests = nrow(z),
        mean_absolute_residual = mean(abs(z$eqc_delta)),
        maximum_absolute_residual = max(abs(z$eqc_delta))
      )
    }))
  } else if (study == "nodes") {
    out <- do.call(rbind, lapply(split(x, x$M), function(z) {
      data.frame(
        M = z$M[1], runs = nrow(z), solver_MAE = mean(abs(z$calibration_residual)),
        fresh_draw_MAE = mean(abs(z$holdout_error)),
        fresh_draw_error_SD = sd(z$holdout_error)
      )
    }))
  } else if (study %in% c("treatment", "treatment-b")) {
    keys <- c("n_items", "target_rho", "te_mean", "te_sd", "te_scale", "method")
    out <- do.call(rbind, lapply(split(x, x[keys], drop = TRUE), function(z) {
      valid <- is.finite(z$p.value)
      n <- sum(valid)
      k <- sum(z$p.value[valid] < .05)
      interval <- if (n) binom.test(k, n)$conf.int else c(NA, NA)
      cbind(z[1, keys, drop = FALSE], data.frame(
        attempts = nrow(z),
        evaluable = n, rejected = k, rejection_rate = if (n) k / n else NA,
        lower = interval[1], upper = interval[2],
        convergence_rate = mean(z$converged)
      ))
    }))
  } else {
    # Studies that already return cell-level summaries retain their full columns.
    out <- x
  }
  write_result(out, repo_path("output", "runs", study, mode, "summary.csv"))
  print(utils::head(out, 12L), row.names = FALSE)
  invisible(out)
}
