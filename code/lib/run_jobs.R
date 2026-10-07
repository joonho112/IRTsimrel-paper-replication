# Resumable, deterministic execution of independent simulation conditions.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT
#
# Each condition writes one checkpoint. A checkpoint is reused only when
# the condition, scientific code, settings, package version, and R version
# match. All study functions set their own seeds; worker order is irrelevant.

bind_rows_union <- function(rows) {
  columns <- unique(unlist(lapply(rows, names)))
  do.call(rbind, lapply(rows, function(x) {
    for (name in setdiff(columns, names(x))) x[[name]] <- NA
    x[columns]
  }))
}
run_jobs <- function(design, fun, context, study, mode, workers = 1L) {
  require_packages(c("digest", "jsonlite"))
  stopifnot(
    nrow(design) > 0L, !anyDuplicated(design$cond_id),
    workers >= 1L, workers == as.integer(workers),
    all(grepl("^[A-Za-z0-9_-]+$", design$cond_id))
  )
  if (.Platform$OS.type == "windows" && workers > 1L) {
    stop("Use --workers 1 on Windows; parallel execution uses forked R processes.")
  }
  out <- repo_path("output", "runs", study, mode)
  dir.create(file.path(out, "checkpoints"), recursive = TRUE, showWarnings = FALSE)
  source_files <- sort(list.files(repo_path("code"),
    pattern = "\\.R$",
    full.names = TRUE, recursive = TRUE
  ))
  source_files <- source_files[grepl("/(lib|01_simulation)/", source_files) | basename(source_files) == "00_setup.R"]
  code_hashes <- vapply(source_files, digest::digest, character(1),
    algo = "sha256", file = TRUE
  )
  version <- as.character(getNamespaceVersion("IRTsimrel"))
  numerical_packages <- c(
    "IRTsimrel", "digest", "jsonlite", "yaml",
    if (study %in% c("recovery", "scores")) "TAM",
    if (study %in% c("treatment", "treatment-b")) "lme4"
  )
  dependencies <- vapply(numerical_packages, function(p) as.character(packageVersion(p)), character(1))
  run_identity <- list(
    study = study, mode = mode, context = context,
    code = unname(code_hashes), package = version,
    dependencies = dependencies, R = as.character(getRversion()), rng = RNGkind()
  )
  do_one <- function(i) {
    row <- design[i, , drop = FALSE]
    signature <- digest::digest(list(run_identity, row), algo = "sha256")
    path <- file.path(out, "checkpoints", paste0(row$cond_id, ".rds"))
    if (file.exists(path)) {
      stored <- readRDS(path)
      if (!identical(stored$signature, signature)) {
        stop(
          "Checkpoint does not match the current settings: ", row$cond_id,
          ". Move the old output directory before starting this configuration."
        )
      }
      return(stored$value)
    }
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    value <- fun(row, context)
    if (!is.data.frame(value) || !nrow(value)) stop("Empty result for ", row$cond_id)
    tmp <- tempfile(tmpdir = dirname(path), fileext = ".rds")
    saveRDS(list(signature = signature, value = value), tmp)
    if (!file.rename(tmp, path)) stop("Could not save checkpoint: ", path)
    value
  }
  chunks <- split(seq_len(nrow(design)), ceiling(seq_len(nrow(design)) / max(1L, workers)))
  values <- vector("list", nrow(design))
  for (indices in chunks) {
    batch <- if (workers == 1L) {
      lapply(indices, do_one)
    } else {
      parallel::mclapply(indices, do_one,
        mc.cores = workers,
        mc.preschedule = FALSE, mc.set.seed = FALSE
      )
    }
    if (any(vapply(batch, inherits, logical(1), "try-error"))) {
      stop("A worker failed; completed checkpoints are retained.")
    }
    values[indices] <- batch
    cat(sprintf("%s: %d / %d conditions complete\n", study, max(indices), nrow(design)))
  }
  result <- bind_rows_union(values)
  write_result(result, file.path(out, "results.csv"))
  write_result(design, file.path(out, "design.csv"))
  jsonlite::write_json(
    list(
      author = "JoonHo Lee (jlee296@ua.edu)",
      study = study, mode = mode, package_version = version, dependencies = as.list(dependencies), conditions = nrow(design),
      rows = nrow(result), workers = workers,
      identity_sha256 = digest::digest(run_identity, algo = "sha256")
    ),
    file.path(out, "run.json"),
    pretty = TRUE, auto_unbox = TRUE
  )
  record_session(file.path(out, "session-info.txt"))
  invisible(result)
}
