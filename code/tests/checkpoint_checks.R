# Verify checkpoint reuse, rejection of changed settings, and safe identifiers.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

source("code/00_setup.R")
use_irtsimrel("0.3.0")
source("code/lib/run_jobs.R")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
run_test <- function() {
  study <- paste0("checkpoint-test-", Sys.getpid())
  out <- repo_path("output", "runs", study)
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  design <- data.frame(cond_id = "condition_1", seed = 101L)
  calls <- 0L
  worker <- function(row, context) {
    calls <<- calls + 1L
    data.frame(cond_id = row$cond_id, answer = context$value)
  }
  a <- run_jobs(design, worker, list(value = 7), study, "smoke")
  b <- run_jobs(design, worker, list(value = 7), study, "smoke")
  stopifnot(identical(a, b), calls == 1L)
  fails <- function(expr) inherits(tryCatch({force(expr); NULL}, error = identity), "error")
  stopifnot(fails(run_jobs(design, worker, list(value = 8), study, "smoke")))
  design$cond_id <- "../escape"
  stopifnot(fails(run_jobs(design, worker, list(value = 7), study, "smoke")))
  stopifnot(fails(use_irtsimrel("0.3.1")), is.null(read_irw_pool(FALSE)))
  cat("Checkpoint and package-isolation checks passed.\n")
}
run_test()
