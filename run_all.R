# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

# Start here from the project directory. Run --help for available commands.
source("code/00_setup.R")
args <- commandArgs(trailingOnly = TRUE)
usage <- function() {
  cat(paste(
    "Rscript run_all.R <command> [--mode smoke|full] [--pool parametric|irw|both] [--workers N]",
    "Commands: verify, quick, example, feasibility, summarize, calibration, algorithms, superpopulation,",
    "          nodes, steps, sampling, recovery, recipes, treatment, treatment-b, dif, scores",
    "Defaults: mode=smoke, pool=parametric, workers=1. See README.md before full runs.",
    sep = "\n"
  ), "\n")
}
if (!length(args) || args[1] %in% c("--help", "-h")) {
  usage()
  quit(status = 0L)
}
command <- args[1]
rest <- args[-1]
if (length(rest) %% 2L) stop("Each option needs one value. See --help.")
opts <- list(mode = "smoke", pool = "parametric", workers = "1", study = "calibration")
if (length(rest)) {
  for (i in seq.int(1L, length(rest), 2L)) {
    key <- sub("^--", "", rest[i])
    if (!startsWith(rest[i], "--") || !key %in% names(opts)) stop("Unknown option: ", rest[i])
    opts[[key]] <- rest[i + 1L]
  }
}
stopifnot(opts$mode %in% c("smoke", "full"), opts$pool %in% c("parametric", "irw", "both"))
workers <- suppressWarnings(as.integer(opts$workers))
if (is.na(workers) || workers < 1L || as.character(workers) != opts$workers) {
  stop("--workers must be a positive integer.")
}
if (command == "summarize") {
  source("code/02_analysis/summarize_run.R")
  summarize_run(opts$study, opts$mode)
} else if (command == "feasibility") {
  source("code/01_simulation/feasibility.R")
  run_feasibility(opts$mode, opts$pool, workers)
} else if (command == "verify") {
  source("code/03_verify.R")
} else if (command == "quick") {
  source("code/03_verify.R")
  source("code/02_analysis/rebuild.R")
} else if (command == "example") {
  source("code/01_simulation/worked_example.R")
} else if (command %in% c("calibration", "algorithms", "superpopulation", "nodes", "steps", "sampling", "recovery")) {
  source("code/01_simulation/validation.R")
  run_validation(command, opts$mode, opts$pool, workers)
} else if (command %in% c("recipes", "treatment", "treatment-b", "dif", "scores")) {
  source("code/01_simulation/extensions.R")
  run_extension(command, opts$mode, opts$pool, workers)
} else {
  stop("Unknown command: ", command)
}
