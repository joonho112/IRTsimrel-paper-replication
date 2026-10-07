# Arguments and diagnostics for fixed-form and item-population calibration.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

# Translate one design row into the package API. The IRW supplies
# difficulties only; discrimination and guessing follow the stated parametric laws.
cell_common_args <- function(cell, pool_data = NULL) {
  model_arg <- switch(cell$model,
    rasch = "rasch",
    `2pl` = "2pl",
    `3pl_g20` = "3pl",
    `3pl_beta` = "3pl",
    stop("unknown model: ", cell$model)
  )

  ip <- list()
  if (cell$pool == "irw") {
    if (is.null(pool_data)) stop("the IRW arm needs pool_data", call. = FALSE)
    ip$difficulty_params <- list(pool = pool_data)
  }
  if (model_arg == "3pl") {
    ip$guessing_params <- if (cell$model == "3pl_g20") {
      list(distribution = "fixed", value = 0.20)
    } else {
      list(distribution = "beta", shape1 = 5, shape2 = 17)
    }
  }

  list(
    n_items = as.integer(cell$n_items),
    model = model_arg,
    latent_shape = cell$latent_shape,
    item_source = if (cell$pool == "irw") "irw" else "parametric",
    item_params = ip
  )
}

cell_sim_args <- function(cell, pool_data = NULL) {
  common <- cell_common_args(cell, pool_data)
  args <- list(
    n_items = common$n_items, model = common$model,
    source = common$item_source
  )
  ip <- common$item_params
  if (!is.null(ip$difficulty_params)) args$difficulty_params <- ip$difficulty_params
  if (!is.null(ip$guessing_params)) args$guessing_params <- ip$guessing_params
  args
}

S2_ITEM_SUPERPOPULATION_AGGREGATION <- "functional_of_mean_information"
S2_ITEM_SUPERPOPULATION_PRIMARY <- "rho_at_expected_information"
S2_ITEM_SUPERPOPULATION_DESCRIPTIVE <- "mean_random_form_reliability"

validate_s2_estimand_contract <- function(cfg) {
  alg <- tryCatch(cfg$s2$algorithms, error = function(e) NULL)
  if (is.null(alg) || is.null(alg$sac_superpop)) {
    stop("S2 must declare algorithms$sac_superpop.", call. = FALSE)
  }
  sp <- alg$sac_superpop
  required <- list(
    schema_version = 2L,
    metric = "info",
    resample_items = TRUE,
    estimand = "item_superpopulation",
    item_aggregation = S2_ITEM_SUPERPOPULATION_AGGREGATION,
    primary_reliability = S2_ITEM_SUPERPOPULATION_PRIMARY,
    descriptive_reliability = S2_ITEM_SUPERPOPULATION_DESCRIPTIVE
  )
  for (nm in names(required)) {
    got <- sp[[nm]]
    want <- required[[nm]]
    if (!identical(got, want)) {
      stop("S2 sac_superpop must declare `", nm, " = ",
        paste0(want, collapse = ", "), "`; found ",
        if (is.null(got)) "NULL" else paste0(got, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  fixed <- intersect(c("sac_info", "sac_cold", "sac_msem"), names(alg))
  bad <- fixed[vapply(fixed, function(nm) {
    spec <- alg[[nm]]
    isTRUE(spec$resample_items) || !identical(spec$estimand, "fixed_form")
  }, logical(1))]
  if (length(bad)) {
    stop("Fixed-form S2 arms changed scope: ", paste(bad, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Calibrate expected information across random forms under package 0.3.1.
# The fixed-form EQC multiplier is an initialization value, not the target estimand.
run_s2_superpopulation_condition <- function(cell, cfg, eqc_c_star = NA_real_,
                                             pool_data = NULL) {
  validate_s2_estimand_contract(cfg)
  spec <- cfg$s2$algorithms$sac_superpop
  rng <- RNGkind()
  seed <- as.integer(cell$seed) + 1L
  only_n <- as.integer(unlist(spec$only_n_items))
  eligible <- as.integer(cell$n_items) %in% only_n
  pre <- spec$preflight_controls
  eva <- spec$evaluation_controls
  burn_in <- floor(as.integer(spec$n_iter) / 2L)

  out <- list(
    cond_id = as.character(cell$cond_id),
    model = as.character(cell$model),
    latent_shape = as.character(cell$latent_shape),
    pool = as.character(cell$pool),
    n_items = as.integer(cell$n_items),
    target_rho = as.numeric(cell$target_rho),
    sac_superpop_eligible = eligible,
    sac_superpop_attempted = eligible,
    sac_superpop_estimand = as.character(spec$estimand),
    sac_superpop_schema_version = as.integer(spec$schema_version),
    sac_superpop_item_aggregation = as.character(spec$item_aggregation),
    sac_superpop_primary_reliability_name =
      as.character(spec$primary_reliability),
    sac_superpop_descriptive_reliability_name =
      as.character(spec$descriptive_reliability),
    sac_superpop_estimating_equation =
      if (eligible) "scaled_mean_information_residual" else NA_character_,
    sac_superpop_package_version =
      as.character(utils::packageVersion("IRTsimrel")),
    sac_superpop_r_version = R.version.string,
    sac_superpop_design_seed = as.integer(cell$seed),
    sac_superpop_seed = seed,
    sac_superpop_master_seed = seed,
    sac_superpop_preflight_stream_seed = NA_integer_,
    sac_superpop_evaluation_stream_seed = NA_integer_,
    sac_superpop_rng_kind = as.character(rng[[1L]]),
    sac_superpop_rng_normal_kind = as.character(rng[[2L]]),
    sac_superpop_rng_sample_kind = as.character(rng[[3L]]),
    sac_superpop_status = if (eligible) "pending" else "not_in_arm",
    sac_superpop_error_class = NA_character_,
    sac_superpop_message = NA_character_,
    sac_superpop_warning_count = 0L,
    sac_superpop_warnings = NA_character_,
    sac_superpop_M_pre = as.integer(spec$M_pre),
    sac_superpop_M_per_iter = as.integer(spec$M_per_iter),
    sac_superpop_n_iter = as.integer(spec$n_iter),
    sac_superpop_burn_in = as.integer(burn_in),
    sac_superpop_calibration_theta_draws_planned =
      as.numeric(spec$M_per_iter) * as.numeric(spec$n_iter),
    sac_superpop_calibration_iterations_completed = 0L,
    sac_superpop_n_post_burn = NA_integer_,
    sac_superpop_preflight_n_forms = as.integer(pre$n_forms),
    sac_superpop_preflight_M = as.integer(pre$M),
    sac_superpop_evaluation_n_forms = as.integer(eva$n_forms),
    sac_superpop_evaluation_M_per_form = as.integer(eva$M),
    sac_superpop_evaluation_theta_draws_planned =
      as.numeric(eva$n_forms) * as.numeric(eva$M),
    sac_superpop_evaluation_theta_draws_completed = 0,
    sac_superpop_n_form_information = 0L,
    sac_superpop_n_finite_form_information = 0L,
    sac_superpop_n_form_reliability = 0L,
    sac_superpop_n_finite_form_reliability = 0L,
    sac_superpop_requested_c_init = if (is.finite(eqc_c_star)) {
      as.numeric(eqc_c_star)
    } else {
      NA_real_
    },
    sac_superpop_warm_start_source = if (is.finite(eqc_c_star)) {
      "archived_eqc_c_star"
    } else {
      "package_default"
    },
    sac_superpop_requested_c_lower = as.numeric(spec$c_bounds[[1L]]),
    sac_superpop_requested_c_upper = as.numeric(spec$c_bounds[[2L]]),
    sac_superpop_effective_c_lower = NA_real_,
    sac_superpop_effective_c_upper = NA_real_,
    sac_superpop_c_star = NA_real_,
    sac_superpop_achieved = NA_real_,
    sac_superpop_delta = NA_real_,
    sac_superpop_achieved_se = NA_real_,
    sac_superpop_projection_rate = NA_real_,
    sac_superpop_projection_side = NA_character_,
    sac_superpop_n_lower_projections = NA_integer_,
    sac_superpop_n_upper_projections = NA_integer_,
    sac_superpop_converged = NA,
    sac_superpop_c_init = NA_real_,
    sac_superpop_init_method = NA_character_,
    sac_superpop_branch = NA_character_,
    sac_superpop_expected_information = NA_real_,
    sac_superpop_rho_at_expected_information = NA_real_,
    sac_superpop_mean_random_form_reliability = NA_real_,
    sac_superpop_jensen_gap = NA_real_,
    sac_superpop_theta_var = NA_real_,
    sac_superpop_elapsed_sec = 0
  )
  if (!eligible) {
    return(as.data.frame(out,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  common <- cell_common_args(cell, pool_data)
  args <- c(common, list(
    target_rho = as.numeric(cell$target_rho),
    reliability_metric = spec$metric,
    M_pre = as.integer(spec$M_pre),
    M_per_iter = as.integer(spec$M_per_iter),
    n_iter = as.integer(spec$n_iter),
    c_bounds = unlist(spec$c_bounds),
    resample_items = TRUE,
    seed = seed,
    verbose = FALSE,
    preflight_controls = pre,
    evaluation_controls = eva
  ), sac_item_aggregation_args(spec))
  if (identical(spec$warm_start_from, "eqc") && is.finite(eqc_c_star)) {
    args$c_init <- as.numeric(eqc_c_star)
  }

  warnings <- character(0)
  err <- NULL
  t0 <- proc.time()[["elapsed"]]
  result <- tryCatch(
    withCallingHandlers(
      do.call(IRTsimrel::sac_calibrate, args),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      err <<- e
      NULL
    }
  )
  out$sac_superpop_elapsed_sec <- proc.time()[["elapsed"]] - t0
  out$sac_superpop_warning_count <- length(warnings)
  if (length(warnings)) {
    out$sac_superpop_warnings <- paste(unique(warnings), collapse = " | ")
  }

  if (is.null(result)) {
    msg <- conditionMessage(err)
    out$sac_superpop_error_class <- paste(class(err), collapse = "|")
    out$sac_superpop_message <- msg
    out$sac_superpop_status <-
      if (grepl("no resolved feasible root", msg, fixed = TRUE)) {
        "unreachable_target"
      } else if (grepl("not attached to one resolved increasing branch",
        msg,
        fixed = TRUE
      )) {
        "no_increasing_branch"
      } else {
        "error"
      }
    return(as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE))
  }

  sp <- extract_item_superpopulation_result(result, spec$item_aggregation)
  dist <- result$achieved_distribution
  prov <- result$rng_provenance
  prov_rng <- as.character(prov$rng_kind)
  out$sac_superpop_schema_version <- as.integer(result$schema_version)
  out$sac_superpop_preflight_stream_seed <- as.integer(prov$preflight_stream_seed)
  out$sac_superpop_evaluation_stream_seed <- as.integer(prov$evaluation_stream_seed)
  out$sac_superpop_rng_kind <- prov_rng[[1L]]
  out$sac_superpop_rng_normal_kind <- prov_rng[[2L]]
  out$sac_superpop_rng_sample_kind <- prov_rng[[3L]]
  out$sac_superpop_status <- as.character(result$calibration_status)
  out$sac_superpop_calibration_iterations_completed <- length(result$trajectory)
  out$sac_superpop_n_post_burn <- as.integer(result$convergence$n_post_burn)
  out$sac_superpop_preflight_n_forms <-
    as.integer(result$calibration_design$preflight_n_forms)
  out$sac_superpop_preflight_M <-
    as.integer(result$calibration_design$preflight_M)
  out$sac_superpop_evaluation_n_forms <- as.integer(dist$n_forms)
  out$sac_superpop_evaluation_M_per_form <- as.integer(dist$M_per_form)
  out$sac_superpop_evaluation_theta_draws_completed <- as.numeric(result$M_final)
  out$sac_superpop_n_form_information <- length(dist$form_mean_information)
  out$sac_superpop_n_finite_form_information <-
    sum(is.finite(dist$form_mean_information))
  out$sac_superpop_n_form_reliability <- length(dist$form_reliability)
  out$sac_superpop_n_finite_form_reliability <-
    sum(is.finite(dist$form_reliability))
  out$sac_superpop_effective_c_lower <- as.numeric(result$c_bounds[[1L]])
  out$sac_superpop_effective_c_upper <- as.numeric(result$c_bounds[[2L]])
  out$sac_superpop_c_star <- as.numeric(result$c_star)
  out$sac_superpop_achieved <- sp$rho_at_expected_information
  out$sac_superpop_delta <- sp$rho_at_expected_information - cell$target_rho
  out$sac_superpop_achieved_se <- as.numeric(result$achieved_se)
  out$sac_superpop_projection_rate <- as.numeric(result$projection_rate)
  lower <- sum(result$projection_side == "lower", na.rm = TRUE)
  upper <- sum(result$projection_side == "upper", na.rm = TRUE)
  out$sac_superpop_n_lower_projections <- as.integer(lower)
  out$sac_superpop_n_upper_projections <- as.integer(upper)
  out$sac_superpop_projection_side <- if (lower && upper) "both" else if (lower) "lower" else if (upper) "upper" else "none"
  out$sac_superpop_converged <- isTRUE(result$convergence$converged)
  out$sac_superpop_c_init <- as.numeric(result$c_init)
  out$sac_superpop_init_method <- as.character(result$init_method)
  out$sac_superpop_branch <- if (isTRUE(result$branch$lost)) {
    "lost"
  } else {
    as.character(result$branch$selected_root_record$direction)
  }
  out$sac_superpop_expected_information <- sp$expected_information
  out$sac_superpop_rho_at_expected_information <-
    sp$rho_at_expected_information
  out$sac_superpop_mean_random_form_reliability <-
    sp$mean_random_form_reliability
  out$sac_superpop_jensen_gap <- sp$jensen_gap
  out$sac_superpop_theta_var <- sp$theta_var
  out$sac_superpop_estimating_equation <- sp$estimating_equation
  as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

sac_item_aggregation_args <- function(spec,
                                      sac_fun = IRTsimrel::sac_calibrate) {
  if (!isTRUE(spec$resample_items)) {
    return(list())
  }
  aggregation <- spec$item_aggregation
  if (!identical(aggregation, S2_ITEM_SUPERPOPULATION_AGGREGATION)) {
    stop("A resampled-items S2 arm must request item_aggregation = \"",
      S2_ITEM_SUPERPOPULATION_AGGREGATION, "\".",
      call. = FALSE
    )
  }
  if (!"item_aggregation" %in% names(formals(sac_fun))) {
    stop("The loaded IRTsimrel package predates the explicit item-aggregation ",
      "contract. Install the reviewed v0.3.1 build before running the S2 ",
      "item-superpopulation arm.",
      call. = FALSE
    )
  }
  list(item_aggregation = aggregation)
}

extract_item_superpopulation_result <- function(
  result,
  declared_aggregation = S2_ITEM_SUPERPOPULATION_AGGREGATION,
  tolerance = 1e-10
) {
  required <- c(
    "schema_version", "item_aggregation", "estimating_equation",
    "achieved_rho", "achieved_mean_information",
    "mean_form_reliability", "theta_var",
    "achieved_distribution"
  )
  absent <- setdiff(required, names(result))
  if (length(absent)) {
    stop("IRTsimrel SAC result lacks the v0.3.1 estimand fields: ",
      paste(absent, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (as.integer(result$schema_version) < 2L ||
    !identical(result$item_aggregation, declared_aggregation) ||
    !identical(
      result$estimating_equation,
      "scaled_mean_information_residual"
    )) {
    stop("IRTsimrel SAC result does not implement the declared ",
      "functional-of-mean-information contract.",
      call. = FALSE
    )
  }
  dist <- result$achieved_distribution
  if (!all(c("form_mean_information", "form_reliability") %in% names(dist))) {
    stop("IRTsimrel SAC achieved_distribution lacks form-level information ",
      "or reliability.",
      call. = FALSE
    )
  }
  two <- item_superpopulation_reliability(
    dist$form_mean_information,
    as.numeric(result$theta_var),
    dist$form_reliability
  )
  identities <- c(
    abs(result$achieved_mean_information - two$expected_information),
    abs(result$achieved_rho - two$rho_at_expected_information),
    abs(result$mean_form_reliability - two$mean_random_form_reliability)
  )
  if (any(!is.finite(identities)) || max(identities) > tolerance) {
    stop("IRTsimrel SAC result fields disagree with its form distribution ",
      "under the declared aggregation.",
      call. = FALSE
    )
  }
  if (two$jensen_gap < -tolerance) {
    stop("Reliability at expected information is below mean random-form ",
      "reliability, contradicting the concavity contract.",
      call. = FALSE
    )
  }
  list(
    expected_information = two$expected_information,
    rho_at_expected_information = two$rho_at_expected_information,
    mean_random_form_reliability = two$mean_random_form_reliability,
    jensen_gap = two$jensen_gap,
    theta_var = as.numeric(result$theta_var),
    item_aggregation = result$item_aggregation,
    estimating_equation = result$estimating_equation
  )
}

# Run the fixed-form solver comparisons for one request. Keep infeasible
# requests and solver failures in the output so success rates retain their denominators.
run_condition <- function(cell, cfg, pool_data = NULL) {
  validate_s2_estimand_contract(cfg)
  alg <- cfg$s2$algorithms
  common <- cell_common_args(cell, pool_data)
  seed0 <- as.integer(cell$seed)

  out <- data.frame(
    cond_id = cell$cond_id, model = cell$model, latent_shape = cell$latent_shape,
    pool = cell$pool, n_items = cell$n_items, target_rho = cell$target_rho,
    stringsAsFactors = FALSE
  )

  timed <- function(expr) {
    t0 <- proc.time()[["elapsed"]]
    val <- tryCatch(force(expr), error = function(e) structure(list(), error = conditionMessage(e)))
    list(value = val, elapsed = proc.time()[["elapsed"]] - t0)
  }
  scalar_chr <- function(v) {
    if (is.null(v) || length(v) == 0L) {
      return(NA_character_)
    }
    as.character(v)[1]
  }
  scalar_num <- function(v) {
    if (is.null(v) || length(v) == 0L) {
      return(NA_real_)
    }
    as.numeric(v)[1]
  }
  got <- function(x, nm, default = NA) {
    v <- tryCatch(x[[nm]], error = function(e) NULL)
    if (is.null(v) || length(v) == 0L) {
      return(default)
    }
    if (length(v) > 1L) {
      return(v[[1]])
    }
    v
  }

  e <- timed(do.call(IRTsimrel::eqc_calibrate, c(common, list(
    target_rho = cell$target_rho, reliability_metric = alg$eqc$metric,
    M = as.integer(alg$eqc$M), c_bounds = unlist(alg$eqc$c_bounds),
    tol = 1e-6, root_policy = alg$eqc$root_policy,
    seed = seed0, verbose = FALSE
  ))))
  ev <- e$value
  eqc_ok <- is.null(attr(ev, "error")) && !is.null(got(ev, "c_star", NULL))
  out$eqc_status <- if (eqc_ok) got(ev$misc, "calibration_status", "ok") else "error"
  out$eqc_message <- if (eqc_ok) NA_character_ else attr(ev, "error")
  out$eqc_c_star <- if (eqc_ok) got(ev, "c_star") else NA_real_
  out$eqc_achieved <- if (eqc_ok) got(ev, "achieved_rho") else NA_real_
  out$eqc_delta <- out$eqc_achieved - cell$target_rho
  out$eqc_root_status <- if (eqc_ok) got(ev$misc, "root_status") else NA_character_
  out$eqc_root_count <- if (eqc_ok) got(ev$misc, "root_count") else NA_integer_
  out$eqc_admissible_roots <- if (eqc_ok) got(ev$misc, "admissible_root_count") else NA_integer_
  out$eqc_topology_status <- if (eqc_ok) got(ev$misc, "topology_status") else NA_character_
  out$eqc_topology_evals <- if (eqc_ok) got(ev$misc, "topology_evaluations") else NA_integer_
  out$eqc_rho_max <- if (eqc_ok) got(ev$misc, "rho_max") else NA_real_
  out$eqc_c_at_rho_max <- if (eqc_ok) got(ev$misc, "c_at_rho_max") else NA_real_
  out$eqc_theta_var <- if (eqc_ok) got(ev, "theta_var") else NA_real_
  out$eqc_elapsed_sec <- e$elapsed

  if (eqc_ok) {
    lam <- got(ev, "lambda_scaled", NULL)
    bet <- got(ev, "beta_vec", NULL)
    gss <- tryCatch(ev$guessing_vec, error = function(e) NULL)
    th <- tryCatch(ev$theta_quad, error = function(e) NULL)
    if (!is.null(lam) && !is.null(bet) && !is.null(th)) {
      lam <- ev$lambda_scaled
      bet <- ev$beta_vec
      th <- as.numeric(ev$theta_quad)
      integrable <- cell$latent_shape != "heavy_tail"

      r_pkg <- rho_three(1, th, bet, lam,
        guessing = gss,
        theta_var = ev$theta_var, tail_integrable = integrable
      )
      out$kernel_vs_package <- r_pkg$rho_tilde - out$eqc_achieved

      r3 <- rho_three(1, th, bet, lam,
        guessing = gss, theta_var = 1,
        tail_integrable = integrable
      )
      out$rho_tilde_at_cstar <- r3$rho_tilde
      out$rho_psd_at_cstar <- r3$rho_psd
      out$w_bar_at_cstar <- r3$w_bar
      out$theta_var_shift <- r3$rho_tilde - r_pkg$rho_tilde
      out$lambda_median <- stats::median(lam)
      out$lambda_p10 <- stats::quantile(lam, .10, names = FALSE)
      out$lambda_p90 <- stats::quantile(lam, .90, names = FALSE)
    }
  }

  sac_arm <- function(spec, prefix) {
    out[[paste0(prefix, "_estimand")]] <<- as.character(spec$estimand)
    if (isTRUE(spec$resample_items)) {
      out[[paste0(prefix, "_item_aggregation")]] <<-
        as.character(spec$item_aggregation)
      out[[paste0(prefix, "_primary_reliability_name")]] <<-
        as.character(spec$primary_reliability)
      out[[paste0(prefix, "_descriptive_reliability_name")]] <<-
        as.character(spec$descriptive_reliability)
    }
    if (!is.null(spec$exclude_shapes) &&
      cell$latent_shape %in% unlist(spec$exclude_shapes)) {
      out[[paste0(prefix, "_status")]] <<- "excluded_by_design"
      return(invisible(NULL))
    }
    if (!is.null(spec$only_n_items) &&
      !(cell$n_items %in% unlist(spec$only_n_items))) {
      out[[paste0(prefix, "_status")]] <<- "not_in_arm"
      return(invisible(NULL))
    }
    aggregation_args <- sac_item_aggregation_args(spec)
    args <- c(common, list(
      target_rho = cell$target_rho, reliability_metric = spec$metric,
      M_pre = as.integer(spec$M_pre), M_per_iter = as.integer(spec$M_per_iter),
      n_iter = as.integer(spec$n_iter), c_bounds = unlist(spec$c_bounds),
      resample_items = isTRUE(spec$resample_items),
      seed = seed0 + 1L, verbose = FALSE
    ), aggregation_args)
    if (!is.null(spec$preflight_controls)) {
      args$preflight_controls <- spec$preflight_controls
    }
    if (!is.null(spec$evaluation_controls)) {
      args$evaluation_controls <- spec$evaluation_controls
    }
    if (identical(spec$warm_start_from, "eqc") && isTRUE(eqc_ok)) {
      args$c_init <- out$eqc_c_star
    }
    s <- timed(do.call(IRTsimrel::sac_calibrate, args))
    sv <- s$value
    ok <- is.null(attr(sv, "error")) && !is.null(got(sv, "c_star", NULL))

    msg <- if (ok) NA_character_ else attr(sv, "error")
    status <- if (ok) {
      got(sv, "calibration_status", "ok")
    } else if (!is.na(msg) && grepl("no resolved feasible root", msg, fixed = TRUE)) {
      "unreachable_target"
    } else if (!is.na(msg) && grepl("not attached to one resolved increasing branch",
      msg,
      fixed = TRUE
    )) {
      "no_increasing_branch"
    } else {
      "error"
    }
    out[[paste0(prefix, "_status")]] <<- status
    out[[paste0(prefix, "_message")]] <<- msg
    out[[paste0(prefix, "_c_star")]] <<- if (ok) scalar_num(got(sv, "c_star", NULL)) else NA_real_
    out[[paste0(prefix, "_achieved")]] <<- if (ok) scalar_num(got(sv, "achieved_rho", NULL)) else NA_real_
    out[[paste0(prefix, "_delta")]] <<- if (ok) {
      scalar_num(got(sv, "achieved_rho", NULL)) - cell$target_rho
    } else {
      NA_real_
    }
    out[[paste0(prefix, "_projection_rate")]] <<- if (ok) scalar_num(got(sv, "projection_rate", NULL)) else NA_real_
    out[[paste0(prefix, "_projection_side")]] <<- if (ok) scalar_chr(got(sv, "projection_side", NULL)) else NA_character_
    out[[paste0(prefix, "_converged")]] <<- if (ok) {
      isTRUE(tryCatch(sv$convergence$converged, error = function(e) NA))
    } else {
      NA_integer_
    }
    out[[paste0(prefix, "_achieved_se")]] <<- if (ok) scalar_num(got(sv, "achieved_se", NULL)) else NA_real_
    out[[paste0(prefix, "_c_init")]] <<- if (ok) scalar_num(got(sv, "c_init", NULL)) else NA_real_
    out[[paste0(prefix, "_init_method")]] <<- if (ok) scalar_chr(got(sv, "init_method", NULL)) else NA_character_

    out[[paste0(prefix, "_branch")]] <<- if (ok) {
      scalar_chr(tryCatch(sv$branch$status, error = function(e) NULL))
    } else {
      NA_character_
    }
    out[[paste0(prefix, "_elapsed_sec")]] <<- s$elapsed
    if (ok && isTRUE(spec$resample_items)) {
      sp <- extract_item_superpopulation_result(sv, spec$item_aggregation)
      out[[paste0(prefix, "_expected_information")]] <<-
        sp$expected_information
      out[[paste0(prefix, "_rho_at_expected_information")]] <<-
        sp$rho_at_expected_information
      out[[paste0(prefix, "_mean_random_form_reliability")]] <<-
        sp$mean_random_form_reliability
      out[[paste0(prefix, "_jensen_gap")]] <<- sp$jensen_gap
      out[[paste0(prefix, "_theta_var")]] <<- sp$theta_var
      out[[paste0(prefix, "_item_aggregation")]] <<- sp$item_aggregation
      out[[paste0(prefix, "_estimating_equation")]] <<-
        sp$estimating_equation
    }
    invisible(NULL)
  }

  sac_arm(alg$sac_info, "sac_info")
  sac_arm(alg$sac_msem, "sac_msem")
  sac_arm(alg$sac_superpop, "sac_superpop")

  if (!is.null(alg$sac_cold)) sac_arm(alg$sac_cold, "sac_cold")

  out$package_version <- as.character(utils::packageVersion("IRTsimrel"))
  out$total_elapsed_sec <- sum(unlist(out[grepl("_elapsed_sec$", names(out))]),
    na.rm = TRUE
  )
  out
}

reclassify_sac_status <- function(cal) {
  for (a in c("sac_info", "sac_msem", "sac_superpop")) {
    sc <- paste0(a, "_status")
    mc <- paste0(a, "_message")
    if (!all(c(sc, mc) %in% names(cal))) next
    bad <- !is.na(cal[[mc]]) & cal[[sc]] == "error"
    cal[[sc]][bad & grepl("no resolved feasible root", cal[[mc]], fixed = TRUE)] <-
      "unreachable_target"
    cal[[sc]][bad & grepl("not attached to one resolved increasing branch",
      cal[[mc]],
      fixed = TRUE
    )] <- "no_increasing_branch"
  }
  cal
}

POOL_NA_MESSAGE <- "must be a finite numeric vector"

# Retry a failed item draw using the original deterministic seed sequence.
sim_item_params_safe <- function(args, seed, max_tries = 20L) {
  for (k in seq_len(max_tries)) {
    s <- if (k == 1L) {
      as.integer(seed)
    } else {
      as.integer((as.numeric(seed) + k * 7919) %% 2000000)
    }
    out <- tryCatch(do.call(IRTsimrel::sim_item_params, c(args, list(seed = s))),
      error = function(e) e
    )
    if (!inherits(out, "error")) {
      attr(out, "n_tries") <- k
      return(out)
    }
    if (!grepl(POOL_NA_MESSAGE, conditionMessage(out), fixed = TRUE)) stop(out)
  }
  stop("item bank could not be drawn in ", max_tries, " attempts at seed ",
    seed, "; the upstream pool defect is firing far more often than ",
    "measured and the cause should be re-examined rather than the limit ",
    "raised",
    call. = FALSE
  )
}
