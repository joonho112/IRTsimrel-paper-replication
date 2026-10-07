# Obtain the exact IRW difficulty pool from its public package source.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT
#
# Check the bundled input, or independently recover it from the pinned
# public IRW package. A different pool is never substituted.

source("code/00_setup.R")
require_packages("digest")
metadata <- utils::read.csv(repo_path("config", "external_inputs.csv"),
  stringsAsFactors = FALSE
)
spec <- metadata[metadata$id == "irw_pool", ]
target <- repo_path("data", "external", "irw_diff_pool.csv")
bundled <- repo_path("data", "irw", "irw_diff_pool.csv")
if (file.exists(bundled) && !"--download" %in% commandArgs(trailingOnly = TRUE)) target <- bundled
if (file.exists(target)) {
  stopifnot(digest::digest(file = target, algo = "sha256") == spec$sha256)
  cat("The existing IRW pool matches the paper's input.\n")
} else {
  lib <- repo_path(".library", "external")
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  local_package <- if ("--download" %in% commandArgs(trailingOnly = TRUE)) {
    ""
  } else {
    tryCatch(find.package("irw"), error = function(e) "")
  }
  env <- new.env(parent = emptyenv())
  if (nzchar(local_package)) {
    utils::data("diff_long", package = "irw", lib.loc = dirname(local_package), envir = env)
  } else {
    # Read only the data object from the pinned source archive. No external
    # package code or installation is needed to recover this input.
    archive <- tempfile(fileext = ".tar.gz")
    extracted <- tempfile(pattern = "irw-data-")
    dir.create(extracted)
    url <- paste0("https://api.github.com/repos/itemresponsewarehouse/Rpkg/tarball/", spec$commit)
    utils::download.file(url, archive, mode = "wb", method = "libcurl")
    members <- utils::untar(archive, list = TRUE)
    member <- members[grepl("/data/diff_long[.]rda$", members)]
    if (length(member) != 1L || grepl("(^|/)\\.\\.(/|$)", member)) {
      stop("Pinned archive does not contain the expected data member.")
    }
    utils::untar(archive, files = member, exdir = extracted)
    load(file.path(extracted, member), envir = env)
    unlink(archive)
    unlink(extracted, recursive = TRUE)
  }
  if (!exists("diff_long", env, inherits = FALSE)) stop("IRW package has no diff_long data.")
  pool <- as.data.frame(env$diff_long)
  tmp <- tempfile(fileext = ".csv")
  utils::write.csv(pool, tmp, row.names = FALSE)
  observed <- digest::digest(file = tmp, algo = "sha256")
  if (!identical(observed, spec$sha256)) {
    stop(
      "The extracted pool differs from the published input (", observed,
      "). Run this downloader with --download to read commit ", spec$commit, " directly."
    )
  }
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(tmp, target, overwrite = FALSE))
  unlink(tmp)
  cat(sprintf(
    "Verified IRW pool: %d item rows, %d instruments.\n",
    nrow(pool), length(unique(pool$dataset))
  ))
}
