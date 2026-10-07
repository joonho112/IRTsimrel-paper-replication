# Design tables and per-condition entry points for the validation studies.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT
#
# Full designs retain the paper's seeds and numerical settings. Smoke runs
# select a small set of conditions; only repeated-sampling and model-recovery
# smoke runs reduce the number of repetitions, as their run records indicate.

select_pool <- function(design, pool) {
  if (!pool %in% c("parametric", "irw", "both")) stop("Invalid pool.")
  if ("pool" %in% names(design) && pool != "both") {
    design <- design[design$pool == pool, , drop = FALSE]
  }
  design
}
calibration_design <- function(mode, pool) {
  d <- select_pool(read_result("results", "s2", "design_final.csv"), pool)
  if (mode == "smoke") {
    d <- d[d$n_items == 10L & d$target_rho == .70, , drop = FALSE]
  }
  d[order(d$cond_id), , drop = FALSE]
}
nodes_design <- function(cfg, mode) {
  s <- cfg$s3a_eqc_M
  d <- expand.grid(
    model = unlist(s$cells$model),
    latent_shape = unlist(s$cells$latent_shape), n_items = unlist(s$cells$n_items),
    M = unlist(s$M), seed_id = seq_len(s$n_seeds), stringsAsFactors = FALSE
  )
  d$pool <- "parametric"
  d$target_rho <- s$target
  d$cond_id <- sprintf(
    "a_%s_%s_%d_M%d_s%02d", d$model,
    substr(d$latent_shape, 1, 4), d$n_items, d$M, d$seed_id
  )
  if (mode == "smoke") {
    d <- d[d$model == "2pl" & d$n_items == 10 & d$M %in% c(500, 20000) &
      d$seed_id == 1, , drop = FALSE]
  }
  d
}
step_design <- function(cfg, mode) {
  s <- cfg$s3b_sac_step
  centre <- list(a = 1, A = 50, gamma = .67, n_iter = 1000L, M_per_iter = 1000L)
  settings <- list(as.data.frame(centre))
  for (name in names(centre)) {
    for (v in setdiff(unlist(s[[name]]), centre[[name]])) {
      row <- centre
      row[[name]] <- v
      settings[[length(settings) + 1L]] <- as.data.frame(row)
    }
  }
  settings <- unique(do.call(rbind, settings))
  settings$setting_id <- sprintf("cfg%02d", seq_len(nrow(settings)))
  d <- merge(expand.grid(
    model = unlist(s$cells$model),
    latent_shape = unlist(s$cells$latent_shape), n_items = unlist(s$cells$n_items),
    seed_id = seq_len(s$n_seeds), stringsAsFactors = FALSE
  ), settings, by = NULL)
  d$pool <- "parametric"
  d$target_rho <- s$target
  d$cond_id <- sprintf(
    "b_%s_%s_%s_s%02d", d$setting_id, d$model,
    substr(d$latent_shape, 1, 4), d$seed_id
  )
  if (mode == "smoke") {
    d <- d[d$model == "2pl" & d$latent_shape == "normal" &
      d$seed_id == 1 & d$setting_id %in% c("cfg01", "cfg02"), , drop = FALSE]
  }
  d
}
recovery_design <- function(cfg, mode) {
  s <- cfg$s4
  d <- expand.grid(
    model = unlist(s$model), target_rho = unlist(s$target),
    latent_shape = unlist(s$latent_shape), n_items = unlist(s$n_items),
    n_persons = unlist(s$n_persons), stringsAsFactors = FALSE
  )
  d$pool <- "parametric"
  d$cond_id <- sprintf(
    "s4_%s_t%02.0f_%s_I%d_N%d", d$model,
    d$target_rho * 100, substr(d$latent_shape, 1, 4), d$n_items, d$n_persons
  )
  d$seed <- 800000L + seq_len(nrow(d)) * 13L
  if (mode == "smoke") {
    d <- d[d$target_rho == .75 & d$latent_shape == "normal" &
      d$n_items == 20 & d$n_persons == 500, , drop = FALSE]
  }
  d
}
calibrate_one <- function(row, cfg, pool) {
  spec <- cfg$s2$algorithms$eqc
  do.call(IRTsimrel::eqc_calibrate, c(cell_common_args(row, pool), list(
    target_rho = row$target_rho, reliability_metric = spec$metric,
    M = as.integer(spec$M), c_bounds = unlist(spec$c_bounds),
    tol = 1e-6, root_policy = spec$root_policy,
    seed = as.integer(row$seed), verbose = FALSE
  )))
}
calibration_worker <- function(row, context) {
  fit <- calibrate_one(row, context$cfg, context$pool_data)
  # Re-evaluate the returned form on exactly the solver's nodes and variance.
  same <- rho_three(1, as.numeric(fit$theta_quad), fit$beta_vec,
    fit$lambda_scaled,
    guessing = fit$guessing_vec,
    theta_var = fit$theta_var, tail_integrable = row$latent_shape != "heavy_tail"
  )
  data.frame(row,
    eqc_c_star = fit$c_star, eqc_achieved = fit$achieved_rho,
    eqc_delta = fit$achieved_rho - row$target_rho,
    eqc_status = fit$misc$calibration_status,
    kernel_vs_package = same$rho_tilde - fit$achieved_rho,
    package_version = as.character(getNamespaceVersion("IRTsimrel")),
    check.names = FALSE
  )
}
algorithm_worker <- function(row, context) {
  # The random-form arm is run separately under 0.3.1; the fixed-form
  # calls here retain the 0.3.0 arguments used by the paper.
  cfg <- context$cfg
  cfg$s2$algorithms$sac_superpop$only_n_items <- -1L
  run_condition(row, cfg, context$pool_data)
}
superpopulation_worker <- function(row, context) {
  run_s2_superpopulation_condition(row, context$cfg,
    eqc_c_star = row$eqc_c_star, pool_data = context$pool_data
  )
}
original_worker <- function(name) {
  force(name)
  function(row, context) {
    fun <- get(name, envir = .GlobalEnv)
    environment(fun) <- list2env(context, parent = .GlobalEnv)
    fun(row)
  }
}
compare_replay <- function(result, study, out) {
  spec <- switch(study,
    calibration = list(
      file = c("s2", "calibration.csv"),
      columns = c("eqc_c_star", "eqc_achieved", "eqc_delta")
    ),
    nodes = list(
      file = c("s3", "m_sensitivity.csv"),
      columns = c("c_star", "calibration_residual", "holdout_error")
    ),
    algorithms = list(
      file = c("s2", "calibration.csv"),
      columns = c("eqc_c_star", "sac_info_c_star", "sac_msem_c_star")
    ),
    steps = list(
      file = c("s3", "step_sensitivity.csv"),
      columns = c("c_star_sac", "c_star_eqc", "c_bias", "achieved")
    ),
    superpopulation = list(
      file = c("s2", "superpopulation_v031.csv"),
      columns = c("sac_superpop_c_star", "sac_superpop_rho_at_expected_information")
    ),
    NULL
  )
  if (is.null(spec)) {
    return(invisible(NULL))
  }
  reference <- read_result("results", spec$file[1], spec$file[2])
  key <- match(result$cond_id, reference$cond_id)
  stopifnot(!anyNA(key))
  comparisons <- do.call(rbind, lapply(spec$columns, function(column) {
    current <- result[[column]]
    saved <- reference[[column]][key]
    finite <- is.finite(current) & is.finite(saved)
    difference <- ifelse(finite, abs(current - saved), NA_real_)
    missing_agrees <- is.na(current) & is.na(saved)
    data.frame(
      cond_id = result$cond_id, quantity = column,
      current = current, retained = saved, absolute_difference = difference,
      tolerance = 1e-8, matches = (finite & difference <= 1e-8) | missing_agrees
    )
  }))
  write_result(comparisons, file.path(out, "comparison.csv"))
  cat(sprintf(
    "Retained-value comparison: %d / %d within tolerance.\n",
    sum(comparisons$matches), nrow(comparisons)
  ))
  # Keep discrepancies visible. A completed rerun is not automatically an
  # exact historical replay on a different numerical platform.
  if (!all(comparisons$matches)) {
    warning("Some results differ from the retained run. Inspect comparison.csv.")
  }
  invisible(comparisons)
}
run_validation <- function(study, mode, pool, workers) {
  require_packages(c("digest", "jsonlite", "yaml"))
  cfg <- yaml::read_yaml(repo_path("config", "validation_design.yml"))
  version <- if (study == "superpopulation") "0.3.1" else "0.3.0"
  use_irtsimrel(version)
  load_scientific_helpers()
  source(repo_path("code", "lib", "study_workers.R"))
  source(repo_path("code", "lib", "run_jobs.R"))
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  context <- list(cfg = cfg)
  if (study %in% c("calibration", "algorithms", "superpopulation", "sampling")) {
    design <- calibration_design(if (study == "calibration") mode else "full", pool)
    if (study == "algorithms" && mode == "smoke") {
      design <- design[design$model == "2pl" & design$latent_shape == "normal" &
        design$n_items == 20 & design$target_rho %in% c(.70, .90), ]
    }
    if (study == "superpopulation") {
      design <- design[design$n_items %in% c(20, 40), ]
      if (mode == "smoke") {
        design <- design[design$model == "2pl" &
          design$latent_shape == "normal" & design$n_items == 20 & design$target_rho == .70, ]
      }
    }
    if (study == "sampling") {
      saved <- read_result("results", "s2", "calibration.csv")
      design$eqc_c_star <- saved$eqc_c_star[match(design$cond_id, saved$cond_id)]
      design <- design[saved$eqc_status[match(design$cond_id, saved$cond_id)] == "ok", ]
      if (mode == "smoke") {
        design <- design[design$model == "2pl" &
          design$latent_shape == "normal" & design$n_items == 20 & design$target_rho == .70, ]
      }
      context$K <- if (mode == "smoke") 10L else as.integer(cfg$s2$replication$K)
      context$NS <- unlist(cfg$s2$replication$n_persons)
    }
    if (study == "superpopulation") {
      saved <- read_result("results", "s2", "calibration.csv")
      design$eqc_c_star <- saved$eqc_c_star[match(design$cond_id, saved$cond_id)]
    }
    context$pool_data <- read_irw_pool(any(design$pool == "irw"))
    fun <- switch(study,
      calibration = calibration_worker,
      algorithms = algorithm_worker,
      superpopulation = superpopulation_worker,
      sampling = original_worker("sampling_worker_original")
    )
  } else if (study == "nodes") {
    design <- nodes_design(cfg, mode)
    context$sa <- cfg$s3a_eqc_M
    fun <- original_worker("node_worker_original")
  } else if (study == "steps") {
    design <- step_design(cfg, mode)
    fun <- original_worker("step_worker_original")
  } else if (study == "recovery") {
    require_packages("TAM")
    design <- recovery_design(cfg, mode)
    context$s4 <- cfg$s4
    if (mode == "smoke") context$s4$n_reps <- 2L
    fun <- original_worker("recovery_worker_original")
  } else {
    stop("Unknown validation study: ", study)
  }
  result <- run_jobs(design, fun, context, study, mode, workers)
  compare_replay(result, study, repo_path("output", "runs", study, mode))
  invisible(result)
}
