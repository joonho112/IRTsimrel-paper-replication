# Install reporting dependencies and the two frozen calibration packages.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT
#
# This is the only setup command that installs packages. Study scripts only
# check their dependencies. Use --packages-only when CRAN dependencies exist.

source("code/00_setup.R")
args <- commandArgs(trailingOnly = TRUE)
packages <- c(
  "digest", "jsonlite", "yaml", "dplyr", "tidyr", "readr",
  "ggplot2", "scales", "patchwork", "TAM", "lme4", "testthat"
)
if (!"--packages-only" %in% args) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) utils::install.packages(missing, repos = "https://cloud.r-project.org")
}
# Install in fresh processes: R cannot load two versions of a namespace at once.
for (version in c("0.3.0", "0.3.1")) {
  lib <- repo_path(".library", version)
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  archive <- repo_path("data", "software", paste0("IRTsimrel_", version, ".tar.gz"))
  utils::install.packages(archive, lib = lib, repos = NULL, type = "source")
  stopifnot(read.dcf(file.path(lib, "IRTsimrel", "DESCRIPTION"), "Version")[[1]] == version)
}
cat("Frozen packages installed in project-local libraries.\n")
