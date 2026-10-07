# Per-condition computations for the published simulation designs.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

# The entry scripts provide the design and settings used by these workers.
# Keep the calibration nodes, fresh evaluation draws, and response draws
# separate: each quantity answers a different question.

# Hold the item-generation settings fixed while varying the integration
# node count. Evaluate the calibrated form on 200,000 independent draws.
node_worker_original <- function(row) {
  common <- cell_common_args(row, NULL)
  seed <- 500000L + as.integer(strtoi(substr(
    digest::digest(row$cond_id, algo = "sha1"), 1,
    6
  ), 16L)) %% 400000L
  t1 <- proc.time()[["elapsed"]]
  e <- tryCatch(do.call(IRTsimrel::eqc_calibrate, c(common, list(
    target_rho = row$target_rho,
    reliability_metric = "info", M = as.integer(row$M), c_bounds = c(0.1, 10), tol = 1e-06,
    seed = seed, verbose = FALSE
  ))), error = function(err) NULL)
  el <- proc.time()[["elapsed"]] - t1
  if (is.null(e)) {
    return(data.frame(
      row[, c(
        "cond_id", "model", "latent_shape", "n_items", "M", "seed_id",
        "target_rho"
      )],
      status = "error", c_star = NA_real_, calibration_residual = NA_real_,
      holdout_error = NA_real_, elapsed_sec = el, stringsAsFactors = FALSE
    ))
  }
  hold <- as.numeric(IRTsimrel::sim_latentG(as.integer(sa$M_holdout),
    shape = row$latent_shape,
    seed = seed + 777777L
  )$theta)
  r <- rho_three(1, hold, e$beta_vec, e$lambda_scaled, guessing = e$guessing_vec, tail_integrable = row$latent_shape !=
    "heavy_tail")
  data.frame(row[, c("cond_id", "model", "latent_shape", "n_items", "M", "seed_id", "target_rho")],
    status = "ok", c_star = e$c_star, calibration_residual = e$achieved_rho - row$target_rho,
    holdout_error = r$rho_tilde - row$target_rho, elapsed_sec = el, stringsAsFactors = FALSE
  )
}

# Change one stochastic-approximation setting at a time around the
# central configuration. Both solvers use the original condition-specific seeds.
step_worker_original <- function(row) {
  common <- cell_common_args(row, NULL)
  seed <- 700000L + as.integer(strtoi(substr(
    digest::digest(row$cond_id, algo = "sha1"), 1,
    6
  ), 16L)) %% 400000L
  e <- tryCatch(do.call(IRTsimrel::eqc_calibrate, c(common, list(
    target_rho = row$target_rho,
    reliability_metric = "info", M = 20000L, c_bounds = c(0.1, 10), tol = 1e-06, seed = seed,
    verbose = FALSE
  ))), error = function(err) NULL)
  t1 <- proc.time()[["elapsed"]]
  s <- tryCatch(do.call(IRTsimrel::sac_calibrate, c(common, list(
    target_rho = row$target_rho,
    reliability_metric = "info", M_pre = 20000L, M_per_iter = as.integer(row$M_per_iter),
    n_iter = as.integer(row$n_iter), c_bounds = c(0.05, 20), resample_items = FALSE, step_params = list(
      a = row$a,
      A = row$A, gamma = row$gamma
    ), seed = seed + 1L, verbose = FALSE
  ))), error = function(err) NULL)
  el <- proc.time()[["elapsed"]] - t1
  data.frame(
    row[, c(
      "cond_id", "setting_id", "model", "latent_shape", "n_items", "seed_id",
      "a", "A", "gamma", "n_iter", "M_per_iter", "target_rho"
    )],
    status = if (is.null(s)) {
      "error"
    } else {
      "ok"
    }, c_star_sac = if (is.null(s)) {
      NA_real_
    } else {
      s$c_star
    }, c_star_eqc = if (is.null(e)) {
      NA_real_
    } else {
      e$c_star
    }, c_bias = if (is.null(s) || is.null(e)) {
      NA_real_
    } else {
      s$c_star - e$c_star
    }, achieved = if (is.null(s)) {
      NA_real_
    } else {
      s$achieved_rho
    }, delta = if (is.null(s)) {
      NA_real_
    } else {
      s$achieved_rho - row$target_rho
    }, projection_rate = if (is.null(s)) {
      NA_real_
    } else {
      (if (length(s$projection_rate)) {
        s$projection_rate[1]
      } else {
        NA_real_
      })
    }, converged = if (is.null(s)) {
      NA
    } else {
      isTRUE(tryCatch(s$convergence$converged, error = function(e) NA))
    }, elapsed_sec = el,
    stringsAsFactors = FALSE
  )
}

# Reconstruct the original calibrated bank before drawing new persons.
# Drawing a bank directly with the same seed would not reproduce the EQC bank.
sampling_worker_original <- function(row) {
  alg <- cfg$s2$algorithms$eqc
  common <- cell_common_args(row, if (row$pool == "irw") {
    pool_data
  } else {
    NULL
  })
  e <- do.call(IRTsimrel::eqc_calibrate, c(common, list(
    target_rho = row$target_rho, reliability_metric = alg$metric,
    M = as.integer(alg$M), c_bounds = unlist(alg$c_bounds), tol = 1e-06, root_policy = alg$root_policy,
    seed = as.integer(row$seed), verbose = FALSE
  )))
  if (!isTRUE(abs(e$c_star - row$eqc_c_star) < 1e-08)) {
    stop("bank reproduction failed for ", row$cond_id, ": c* ", signif(e$c_star, 8), " against the stored ",
      signif(row$eqc_c_star, 8),
      call. = FALSE
    )
  }
  beta <- e$beta_vec
  lam <- e$lambda_scaled
  g <- e$guessing_vec
  integrable <- row$latent_shape != "heavy_tail"
  do.call(rbind, lapply(NS, function(n) {
    vals <- vapply(seq_len(K), function(k) {
      th <- as.numeric(IRTsimrel::sim_latentG(n, shape = row$latent_shape, seed = as.integer((row$seed +
        k * 7L + n) %% 2000000L))$theta)
      r <- rho_three(1, th, beta, lam, guessing = g, weights = NULL, tail_integrable = integrable)
      c(r$rho_tilde, r$rho_psd, if (is.na(r$w_bar)) r$w_bar_truncated else r$w_bar)
    }, numeric(3))
    data.frame(
      cond_id = row$cond_id, model = row$model, latent_shape = row$latent_shape,
      pool = row$pool, n_items = row$n_items, target_rho = row$target_rho, n_persons = n,
      K = K, tilde_mean = mean(vals[1, ]), tilde_sd = stats::sd(vals[1, ]), tilde_q025 = stats::quantile(vals[1, ], 0.025, names = FALSE), tilde_q975 = stats::quantile(vals[1, ], 0.975, names = FALSE),
      psd_mean = mean(vals[2, ]), psd_sd = stats::sd(vals[2, ]), wbar_mean = mean(vals[3, ]), wbar_sd = stats::sd(vals[3, ]), wbar_defined = integrable, population_tilde = e$achieved_rho,
      stringsAsFactors = FALSE
    )
  }))
}

# Draw binary responses from the calibrated form and fit the specified
# TAM model. Analytic information and fitted score reliability answer different questions.
recovery_worker_original <- function(row) {
  common <- cell_common_args(row, NULL)
  e <- tryCatch(do.call(IRTsimrel::eqc_calibrate, c(common, list(
    target_rho = row$target_rho,
    reliability_metric = "info", M = 20000L, c_bounds = c(0.1, 10), tol = 1e-06, seed = as.integer(row$seed),
    verbose = FALSE
  ))), error = function(err) NULL)
  if (is.null(e)) {
    return(data.frame(row[, c(
      "cond_id", "model", "target_rho", "latent_shape", "n_items",
      "n_persons"
    )], status = "calibration_error", stringsAsFactors = FALSE))
  }
  beta0 <- e$beta_vec
  lam0 <- e$lambda_scaled
  analytic <- rho_three(1, as.numeric(e$theta_quad), beta0, lam0,
    guessing = e$guessing_vec,
    theta_var = 1
  )
  reps <- lapply(seq_len(s4$n_reps), function(k) {
    seed_k <- as.integer((row$seed + k * 31L) %% 2000000L)
    th <- as.numeric(IRTsimrel::sim_latentG(row$n_persons, shape = row$latent_shape, seed = seed_k)$theta)
    p <- stats::plogis(outer(th, beta0, "-") * rep(lam0, each = length(th)))
    set.seed(seed_k + 1L)
    y <- matrix(stats::rbinom(length(p), 1L, p), nrow = length(th))
    keep <- apply(y, 2, function(col) length(unique(col)) > 1L)
    if (sum(keep) < 3L) {
      return(NULL)
    }
    fit <- tryCatch(if (row$model == "rasch") {
      TAM::tam.mml(y[, keep, drop = FALSE], irtmodel = "1PL", control = list(progress = FALSE))
    } else {
      TAM::tam.mml.2pl(y[, keep, drop = FALSE], control = list(progress = FALSE))
    }, error = function(err) NULL)
    if (is.null(fit)) {
      return(NULL)
    }
    eap <- tryCatch(as.numeric(fit$EAP.rel), error = function(err) NA_real_)
    wle <- tryCatch(as.numeric(TAM::tam.wle(fit, progress = FALSE)$WLE.rel[1]), error = function(err) NA_real_)
    lam_hat <- if (row$model == "rasch") {
      rep(1, sum(keep))
    } else {
      tryCatch(as.numeric(fit$B[, 2, 1]), error = function(err) rep(NA_real_, sum(keep)))
    }
    b_hat <- tryCatch(as.numeric(fit$xsi$xsi), error = function(err) rep(NA_real_, sum(keep)))
    lam_true <- lam0[keep]
    b_true <- beta0[keep]
    align <- if (row$model == "rasch") {
      1
    } else {
      exp(mean(log(lam_true[lam_hat > 0]), na.rm = TRUE) - mean(log(lam_hat[lam_hat >
        0]), na.rm = TRUE))
    }
    c(
      lam_cor = suppressWarnings(stats::cor(lam_hat * align, lam_true, use = "complete.obs")),
      lam_rmse = sqrt(mean((lam_hat * align - lam_true)^2, na.rm = TRUE)), b_cor = suppressWarnings(stats::cor(b_hat,
        b_true,
        use = "complete.obs"
      )), b_rmse = sqrt(mean((b_hat - b_true)^2, na.rm = TRUE)),
      eap_rel = eap, wle_rel = wle
    )
  })
  reps <- reps[!vapply(reps, is.null, logical(1))]
  if (!length(reps)) {
    return(data.frame(row[, c(
      "cond_id", "model", "target_rho", "latent_shape", "n_items",
      "n_persons"
    )], status = "all_fits_failed", stringsAsFactors = FALSE))
  }
  m <- do.call(rbind, reps)
  data.frame(row[, c("cond_id", "model", "target_rho", "latent_shape", "n_items", "n_persons")],
    status = "ok", n_fits = nrow(m), c_star = e$c_star, rho_tilde = analytic$rho_tilde, rho_psd = analytic$rho_psd,
    w_bar = analytic$w_bar, lam_cor = mean(m[, "lam_cor"], na.rm = TRUE), lam_rmse = mean(m[
      ,
      "lam_rmse"
    ], na.rm = TRUE), b_cor = mean(m[, "b_cor"], na.rm = TRUE), b_rmse = mean(m[
      ,
      "b_rmse"
    ], na.rm = TRUE), eap_rel_mean = mean(m[, "eap_rel"], na.rm = TRUE), eap_rel_sd = stats::sd(m[
      ,
      "eap_rel"
    ], na.rm = TRUE), wle_rel_mean = mean(m[, "wle_rel"], na.rm = TRUE), wle_rel_sd = stats::sd(m[
      ,
      "wle_rel"
    ], na.rm = TRUE), eap_minus_psd = mean(m[, "eap_rel"], na.rm = TRUE) - analytic$rho_psd,
    eap_minus_tilde = mean(m[, "eap_rel"], na.rm = TRUE) - analytic$rho_tilde, wle_minus_wbar = mean(m[
      ,
      "wle_rel"
    ], na.rm = TRUE) - analytic$w_bar, wle_minus_psd = mean(m[, "wle_rel"],
      na.rm = TRUE
    ) - analytic$rho_psd, stringsAsFactors = FALSE
  )
}

# Compare sum scores, EAP, and WLE under balanced and unequal latent spreads.
scores_worker_original <- function(row) {
  ip <- IRTsimrel::sim_item_params(N_ITEMS, model = "2pl", source = "parametric", seed = as.integer(row$seed))
  d <- ip$data
  quad <- stats::rnorm(20000L)
  f <- function(cc) rho_three(cc, quad, d$beta, d$lambda, theta_var = 1)$rho_tilde - TARGET_REF
  cstar <- stats::uniroot(f, c(0.05, 10), tol = 1e-08)$root
  lam <- d$lambda * cstar
  rho_ref <- rho_three(1, stats::rnorm(20000L, 0, 1), d$beta, lam, theta_var = NULL)$rho_tilde
  rho_foc <- rho_three(1, stats::rnorm(20000L, row$delta, row$focal_sd), d$beta, lam, theta_var = NULL)$rho_tilde
  reps <- lapply(seq_len(N_REPS), function(k) {
    seed_k <- as.integer((row$seed + k * 41L) %% 2000000L)
    set.seed(seed_k)
    th_r <- stats::rnorm(N_PER_GROUP, 0, 1)
    th_f <- stats::rnorm(N_PER_GROUP, row$delta, row$focal_sd)
    th <- c(th_r, th_f)
    grp <- rep(c(0L, 1L), each = N_PER_GROUP)
    p <- irt_prob(th, d$beta, lam)
    y <- matrix(stats::rbinom(length(p), 1L, p), nrow = length(th))
    keep <- apply(y, 2, function(col) length(unique(col)) > 1L)
    if (sum(keep) < 5L) {
      return(NULL)
    }
    y <- y[, keep, drop = FALSE]
    sum_score <- rowSums(y)
    fit <- tryCatch(TAM::tam.mml.2pl(y, control = list(progress = FALSE)), error = function(e) NULL)
    if (is.null(fit)) {
      return(NULL)
    }
    eap <- tryCatch(as.numeric(fit$person$EAP), error = function(e) rep(NA_real_, length(th)))
    wle <- tryCatch(as.numeric(TAM::tam.wle(fit, progress = FALSE)$theta), error = function(e) {
      rep(
        NA_real_,
        length(th)
      )
    })
    one <- function(sc, nm) {
      if (all(is.na(sc))) {
        return(NULL)
      }
      tt <- stats::t.test(sc[grp == 1L], sc[grp == 0L])
      sd_ref <- stats::sd(sc[grp == 0L])
      data.frame(score = nm, p = tt$p.value, d_hat = (mean(sc[grp == 1L]) - mean(sc[grp ==
        0L])) / sd_ref, stringsAsFactors = FALSE)
    }
    do.call(rbind, list(one(sum_score, "sum"), one(eap, "EAP"), one(wle, "WLE")))
  })
  reps <- reps[!vapply(reps, is.null, logical(1))]
  if (!length(reps)) {
    return(NULL)
  }
  all <- do.call(rbind, reps)
  do.call(rbind, lapply(split(all, all$score), function(g) {
    data.frame(
      cond_id = row$cond_id,
      delta = row$delta, condition = row$condition, focal_sd = row$focal_sd, score = g$score[1],
      n_reps = nrow(g), rho_reference = rho_ref, rho_focal = rho_foc, rho_gap = rho_ref - rho_foc,
      reject_rate = mean(g$p < ALPHA), mcse = sqrt(mean(g$p < ALPHA) * (1 - mean(g$p < ALPHA)) / nrow(g)),
      mean_d_hat = mean(g$d_hat, na.rm = TRUE), bias_d_hat = mean(g$d_hat, na.rm = TRUE) -
        row$delta, sd_d_hat = stats::sd(g$d_hat, na.rm = TRUE), stringsAsFactors = FALSE
    )
  }))
}

# Keep a numerator and denominator for each independent item form.
# The same named person/response streams pair the variants and calibration arms.
dif_worker_original <- function(row) {
  sid <- row$structural_cell_id_full
  qcal_seed <- s5_seed_lookup(seed_registry, sid, NA_integer_, "quadrature_calibration", NA_integer_)
  qeval_seed <- s5_seed_lookup(seed_registry, sid, NA_integer_, "quadrature_evaluation", NA_integer_)
  quad <- as.numeric(IRTsimrel::sim_latentG(M_CALIB, shape = row$latent_shape, seed = as.integer(qcal_seed$seed))$theta)
  quad_eval <- as.numeric(IRTsimrel::sim_latentG(M_EVAL, shape = row$latent_shape, seed = as.integer(qeval_seed$seed))$theta)
  per_rep <- vector("list", R_REP)
  for (k in seq_len(R_REP)) {
    bank_draw <- s5_draw_item_bank(cell_sim_args(row, if (row$pool == "irw") {
      pool_data
    } else {
      NULL
    }), seed_registry, sid, k, as.integer(cfg$rng$max_bank_attempts))
    d <- bank_draw$value$data
    bank <- list(beta = d$beta, lambda_base = d$lambda, guessing = if (grepl("^3pl", row$model)) d$guessing else NULL)
    cvals <- c(uncontrolled = 1)
    for (arm in ARMS) {
      spec <- base_cfg$arms[[arm]]
      if (identical(spec$kind, "uncontrolled")) {
        next
      }
      target <- as.numeric(spec$target_rho)
      f <- function(cc) {
        rho_three(cc, quad, bank$beta, bank$lambda_base,
          guessing = bank$guessing,
          theta_var = 1
        )$rho_tilde - target
      }
      cvals[arm] <- if (f(ROOT_INTERVAL[1]) > 0 || f(ROOT_INTERVAL[2]) < 0) {
        NA_real_
      } else {
        stats::uniroot(f, ROOT_INTERVAL, tol = ROOT_TOL)$root
      }
    }
    realized <- vapply(ARMS, function(arm) {
      cc <- cvals[[arm]]
      if (!is.finite(cc)) {
        NA_real_
      } else {
        rho_three(cc, quad_eval, bank$beta, bank$lambda_base,
          guessing = bank$guessing,
          theta_var = 1
        )$rho_tilde
      }
    }, numeric(1))
    names(realized) <- ARMS
    usable <- cvals[is.finite(cvals)]
    n_dif <- max(1L, round(nrow(d) * base_cfg$dif$share_of_items))
    dif_items <- round(seq(1, nrow(d), length.out = n_dif + 2L))[-c(1, n_dif + 2L)]
    dif_items <- unique(order(d$beta)[dif_items])
    ability_seed <- s5_seed_lookup(seed_registry, sid, k, "latent_abilities", NA_integer_)
    uniform_seed <- s5_seed_lookup(seed_registry, sid, k, "response_uniforms", NA_integer_)
    form_cluster_id <- s5_sha256_text(s5_canonical_key(list(structural_cell_id = sid, replication = k)))
    vv <- vector("list", nrow(VARIANTS))
    for (j in seq_len(nrow(VARIANTS))) {
      vr <- VARIANTS[j, , drop = FALSE]
      condition_row <- catalog[catalog$structural_cell_id_full == sid & catalog$variant ==
        vr$variant, , drop = FALSE]
      stopifnot(nrow(condition_row) == 1L)
      res <- dif_replicate(bank, usable,
        dif_items = dif_items, dif_shift = base_cfg$dif$beta_shift,
        n_per_group = as.integer(row$n_per_group), seed = as.integer(ability_seed$seed),
        uniform_seed = as.integer(uniform_seed$seed), tests = "mh", dif_scale = vr$dif_scale,
        impact = vr$impact_mean_focal, purify = vr$purify, alpha = as.numeric(cfg$alpha)
      )
      variant_form_id <- s5_sha256_text(s5_canonical_key(list(
        condition_id = condition_row$condition_id_full,
        replication = k
      )))
      vv[[j]] <- s5_aggregate_form(
        res, row, condition_row, base_cfg, cfg, k, form_cluster_id,
        variant_form_id, dif_items, cvals, realized, bank_draw, ability_seed, uniform_seed,
        rng_kind
      )
    }
    per_rep[[k]] <- do.call(rbind, vv)
  }
  do.call(rbind, per_rep)
}
