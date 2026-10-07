# Screen the target ladder and reconstruct the calibration design.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

run_feasibility <- function(mode, pool, workers) {
  use_irtsimrel("0.3.0")
  load_scientific_helpers()
  source(repo_path("code", "lib", "run_jobs.R"))
  source(repo_path("code", "01_simulation", "validation.R"))
  require_packages(c("yaml", "digest"))
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  cfg <- yaml::read_yaml(repo_path("config", "validation_design.yml"))
  s <- cfg$s2
  cells <- expand.grid(
    model = unlist(s$factors$model),
    latent_shape = unlist(s$factors$latent_shape), pool = unlist(s$factors$pool),
    n_items = unlist(s$factors$n_items), stringsAsFactors = FALSE
  )
  cells$cond_id <- vapply(seq_len(nrow(cells)), function(i) {
    substr(digest::digest(paste(cells$model[i], cells$latent_shape[i],
      cells$pool[i], cells$n_items[i],
      sep = "|"
    ), algo = "sha1"), 1L, 10L)
  }, character(1))
  cells <- select_pool(cells, pool)
  if (mode == "smoke") cells <- cells[cells$model == "2pl" & cells$n_items == 20, ]
  context <- list(s = s, pool_data = read_irw_pool(any(cells$pool == "irw")))
  worker <- function(row, ctx) {
    common <- cell_common_args(row, ctx$pool_data)
    fit <- IRTsimrel::check_feasibility(
      n_items = common$n_items, model = common$model,
      latent_shape = common$latent_shape, item_source = common$item_source,
      item_params = common$item_params, c_bounds = unlist(ctx$s$feasibility$c_bounds),
      M = as.integer(ctx$s$feasibility$M),
      seed = as.integer(strtoi(substr(row$cond_id, 1, 6), 16L)) %% 90000L + 1000L,
      verbose = FALSE
    )
    data.frame(row,
      rho_min = fit$rho_range_info[1], rho_max = fit$rho_range_info[2],
      rho_msem_min = fit$rho_range_msem[1], rho_msem_max = fit$rho_range_msem[2]
    )
  }
  result <- run_jobs(cells, worker, context, "feasibility", mode, workers)
  design <- do.call(rbind, lapply(seq_len(nrow(result)), function(i) {
    f <- result[i, ]
    ladder <- unlist(s$target_ladder)
    keep <- ladder[ladder > f$rho_min & ladder < f$rho_max]
    if (!length(keep)) {
      return(NULL)
    }
    data.frame(f[rep(1, length(keep)), c("cond_id", "model", "latent_shape", "pool", "n_items")],
      target_rho = keep, cell_rho_min = f$rho_min, cell_rho_max = f$rho_max, row.names = NULL
    )
  }))
  design$cond_id <- sprintf("%s_t%02.0f", design$cond_id, design$target_rho * 100)
  design$seed <- as.integer(cfg$global_seed %% 90000L + vapply(
    design$cond_id,
    function(x) strtoi(substr(digest::digest(x, algo = "sha1"), 1, 6), 16L) %% 800000L,
    numeric(1)
  )) + 1000L
  write_result(design, repo_path("output", "runs", "feasibility", mode, "accepted_design.csv"))
  saved <- select_pool(read_result("results", "s2", "design_final.csv"), pool)
  if (mode == "smoke") saved <- saved[saved$model == "2pl" & saved$n_items == 20, ]
  stopifnot(
    setequal(design$cond_id, saved$cond_id),
    all(design$seed == saved$seed[match(design$cond_id, saved$cond_id)])
  )
  message("Screen and calibration seeds match the retained design: ", nrow(design), " requests.")
  invisible(design)
}
