# Project paths, dependency checks, and isolated package versions.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT
#
# Source this file before a study. Nothing is downloaded automatically.
# The 0.3.0 and 0.3.1 packages live in different project libraries because
# the paper's fixed-form results and item-population results used different versions.

repo_root <- function(start = getwd()) {
  preset <- getOption("irtsimrel.replication.root")
  if (!is.null(preset)) {
    return(preset)
  }
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "IRTsimrel-paper-replication.Rproj"))) {
      options(irtsimrel.replication.root = path)
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) stop("Run from the replication project directory.", call. = FALSE)
    path <- parent
  }
}
repo_path <- function(...) file.path(repo_root(), ...)
read_result <- function(...) {
  utils::read.csv(repo_path("data", "precomputed", ...),
    stringsAsFactors = FALSE, check.names = FALSE
  )
}
write_result <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
  invisible(path)
}
require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Install the required packages first: ", paste(missing, collapse = ", "),
      ". See code/install_dependencies.R.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
use_irtsimrel <- function(version = "0.3.0", install = FALSE) {
  if (!version %in% c("0.3.0", "0.3.1")) stop("Unsupported package version.")
  lib <- repo_path(".library", version)
  archive <- repo_path("data", "software", paste0("IRTsimrel_", version, ".tar.gz"))
  description <- file.path(lib, "IRTsimrel", "DESCRIPTION")
  if (!file.exists(description)) {
    if (!install) stop("Install the frozen packages first: Rscript code/install_dependencies.R --packages-only", call. = FALSE)
    dir.create(lib, recursive = TRUE, showWarnings = FALSE)
    utils::install.packages(archive, repos = NULL, type = "source", lib = lib)
  }
  installed <- read.dcf(description, fields = "Version")[[1]]
  if (!identical(installed, version)) stop("The project library has the wrong version.")
  if ("IRTsimrel" %in% loadedNamespaces()) {
    loaded <- as.character(getNamespaceVersion("IRTsimrel"))
    if (!identical(loaded, version)) {
      stop(
        "Start a fresh R process before changing IRTsimrel versions.",
        call. = FALSE
      )
    }
    if (!identical(
      normalizePath(getNamespaceInfo(asNamespace("IRTsimrel"), "path")),
      normalizePath(file.path(lib, "IRTsimrel"))
    )) {
      stop("IRTsimrel was loaded from a different library. Start a fresh R process.")
    }
  }
  .libPaths(c(lib, .libPaths()))
  loadNamespace("IRTsimrel", lib.loc = lib)
  invisible(lib)
}
load_scientific_helpers <- function() {
  for (name in c("rng", "reliability", "theory", "calibrate", "implicit", "ilhte", "dif", "dif_cluster")) {
    sys.source(repo_path("code", "lib", paste0("fct_", name, ".R")), envir = .GlobalEnv)
  }
}
read_irw_pool <- function(required = FALSE) {
  if (!required) {
    return(NULL)
  }
  path <- repo_path("data", "irw", "irw_diff_pool.csv")
  if (!file.exists(path)) path <- repo_path("data", "external", "irw_diff_pool.csv")
  if (!file.exists(path)) {
    if (required) {
      stop(
        "The IRW pool is required for these cells. Run Rscript data/fetch_irw_pool.R first.",
        call. = FALSE
      )
    }
    return(NULL)
  }
  require_packages("digest")
  expected <- utils::read.csv(repo_path("config", "external_inputs.csv"),
    stringsAsFactors = FALSE
  )
  wanted <- expected$sha256[expected$id == "irw_pool"]
  observed <- digest::digest(file = path, algo = "sha256")
  if (!identical(observed, wanted)) stop("IRW pool checksum differs from the paper's frozen input.")
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
record_session <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(c(
    "Author: JoonHo Lee (jlee296@ua.edu)",
    capture.output(sessionInfo()),
    paste("RNG:", paste(RNGkind(), collapse = ", "))
  ), path)
}
if (getRversion() < "4.3.0") stop("R 4.3.0 or newer is required.")
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"
)
options(stringsAsFactors = FALSE, scipen = 7)
