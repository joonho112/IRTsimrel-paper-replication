# Form-level denominators, named random seeds, and paired DIF uncertainty.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

s5_sha256_text <- function(x) {
  digest::digest(x, algo = "sha256", serialize = FALSE)
}

s5_canon_value <- function(x) {
  if (length(x) > 1L) {
    return(paste(vapply(x, s5_canon_value, character(1)),
      collapse = ","
    ))
  }
  if (is.logical(x)) {
    return(if (isTRUE(x)) "true" else "false")
  }
  if (is.numeric(x)) {
    return(format(x,
      digits = 17L, scientific = FALSE,
      trim = TRUE
    ))
  }
  utils::URLencode(as.character(x), reserved = TRUE)
}

s5_canonical_key <- function(fields) {
  stopifnot(
    length(fields) > 0L, !is.null(names(fields)),
    all(nzchar(names(fields)))
  )
  paste(sprintf(
    "%s=%s", names(fields),
    vapply(fields, s5_canon_value, character(1))
  ), collapse = "|")
}

s5_variant_table <- function(rerun_cfg) {
  nm <- names(rerun_cfg$variants)
  stopifnot(length(nm) > 0L, !anyDuplicated(nm))
  do.call(rbind, lapply(nm, function(v) {
    x <- rerun_cfg$variants[[v]]
    data.frame(
      variant = v, dif_scale = as.character(x$dif_scale),
      purify = isTRUE(x$purify),
      impact_mean_focal = as.numeric(x$impact_mean_focal),
      stringsAsFactors = FALSE
    )
  }))
}

s5_structural_cells <- function(base_cfg) {
  x <- expand.grid(
    model = unlist(base_cfg$factors$model),
    pool = unlist(base_cfg$factors$pool),
    latent_shape = unlist(base_cfg$factors$latent_shape),
    n_items = as.integer(unlist(base_cfg$factors$n_items)),
    n_per_group =
      as.integer(unlist(base_cfg$factors$n_per_group)),
    stringsAsFactors = FALSE
  )
  x$structural_key <- vapply(seq_len(nrow(x)), function(i) {
    s5_canonical_key(list(
      study = "S5", schema = 1L,
      model = x$model[i], pool = x$pool[i],
      latent_shape = x$latent_shape[i], n_items = x$n_items[i],
      n_per_group = x$n_per_group[i]
    ))
  }, character(1))
  x$structural_cell_id_full <- vapply(
    x$structural_key, s5_sha256_text,
    character(1)
  )

  x$cond_id <- x$structural_cell_id_full
  x$cost <- x$n_items * x$n_per_group
  stopifnot(!anyDuplicated(x$structural_cell_id_full), nrow(x) == 32L)
  x
}

s5_effect_contract <- function(model, variant, dif_scale) {
  is_2pl <- identical(as.character(model), "2pl")
  fixed_lp <- identical(as.character(dif_scale), "logodds")
  if (fixed_lp && is_2pl) {
    return(list(
      effect_definition_id = "2pl_fixed_conditional_logodds_and_lp_0p5",
      effect_scale_contract = "conditional_log_odds_and_linear_predictor",
      effect_interpretation = paste(
        "focal difficulty shift=0.5/(lambda_i*c);",
        "conditional log-odds and linear-predictor shifts are both +0.5"
      )
    ))
  }
  if (fixed_lp && !is_2pl) {
    return(list(
      effect_definition_id = "3pl_fixed_linear_predictor_0p5",
      effect_scale_contract = "linear_predictor_only",
      effect_interpretation = paste(
        "focal difficulty shift=0.5/(lambda_i*c); linear-predictor shift",
        "is +0.5; guessing makes observed conditional log odds nonconstant"
      )
    ))
  }
  if (is_2pl) {
    return(list(
      effect_definition_id = "2pl_fixed_theta_shift_0p5",
      effect_scale_contract = "conditional_log_odds_and_linear_predictor",
      effect_interpretation = paste(
        "focal difficulty shift=+0.5 theta units; conditional log-odds",
        "and linear-predictor shifts equal lambda_i*c*0.5"
      )
    ))
  }
  list(
    effect_definition_id = "3pl_fixed_theta_shift_0p5",
    effect_scale_contract = "linear_predictor_only",
    effect_interpretation = paste(
      "focal difficulty shift=+0.5 theta units; linear-predictor shift",
      "equals lambda_i*c*0.5; no constant observed conditional log-odds claim"
    )
  )
}

s5_condition_catalog <- function(base_cfg, rerun_cfg, structural = NULL) {
  if (is.null(structural)) structural <- s5_structural_cells(base_cfg)
  variants <- s5_variant_table(rerun_cfg)
  out <- merge(structural, variants, by = NULL, sort = FALSE)
  arms <- as.character(unlist(rerun_cfg$arms))
  reps <- as.integer(rerun_cfg$replications_per_structural_cell)
  out$effect_definition_id <- character(nrow(out))
  out$effect_scale_contract <- character(nrow(out))
  out$effect_interpretation <- character(nrow(out))
  out$condition_key <- character(nrow(out))
  out$condition_id_full <- character(nrow(out))
  out$source_cond_id <- character(nrow(out))
  for (i in seq_len(nrow(out))) {
    ec <- s5_effect_contract(out$model[i], out$variant[i], out$dif_scale[i])
    out$effect_definition_id[i] <- ec$effect_definition_id
    out$effect_scale_contract[i] <- ec$effect_scale_contract
    out$effect_interpretation[i] <- ec$effect_interpretation
    out$condition_key[i] <- s5_canonical_key(list(
      study = "S5_form_cluster", schema = rerun_cfg$schema_version,
      structural_cell_id = out$structural_cell_id_full[i],
      variant = out$variant[i], dif_scale = out$dif_scale[i],
      purify = out$purify[i], impact_mean_focal = out$impact_mean_focal[i],
      effect_definition = out$effect_definition_id[i],
      replications = reps, arms = arms,
      conditions = unlist(base_cfg$dif$conditions),
      share_of_items = base_cfg$dif$share_of_items,
      beta_shift = base_cfg$dif$beta_shift,
      tests = unlist(rerun_cfg$tests), alpha = rerun_cfg$alpha,
      matching = rerun_cfg$matching
    ))
    out$condition_id_full[i] <- s5_sha256_text(out$condition_key[i])
    legacy <- if (out$variant[i] == "baseline") {
      paste(out$model[i], out$pool[i], out$latent_shape[i], out$n_items[i],
        out$n_per_group[i],
        sep = "|"
      )
    } else {
      paste(out$variant[i], out$model[i], out$pool[i],
        out$latent_shape[i], out$n_items[i], out$n_per_group[i],
        sep = "|"
      )
    }
    out$source_cond_id[i] <- substr(
      digest::digest(legacy, algo = "sha1"),
      1L, 10L
    )
  }
  stopifnot(!anyDuplicated(out$condition_id_full), nrow(out) == 128L)
  out
}

s5_seed_token <- function(namespace, global_label, structural_cell_id_full,
                          replication = NA_integer_, stream,
                          attempt = NA_integer_) {
  s5_canonical_key(list(
    namespace = namespace, global_label = global_label,
    structural_cell_id = structural_cell_id_full,
    replication = if (is.na(replication)) "condition" else replication,
    stream = stream,
    attempt = if (is.na(attempt)) "none" else attempt
  ))
}

s5_allocate_content_hash_seeds <- function(registry) {
  stopifnot(
    is.data.frame(registry), "seed_token" %in% names(registry),
    !anyDuplicated(registry$seed_token)
  )
  registry$seed_token_sha256 <- vapply(
    registry$seed_token, s5_sha256_text,
    character(1)
  )
  raw <- vapply(registry$seed_token, digest::digest2int, integer(1))
  registry$raw_seed <- as.integer(
    (as.double(raw) + 2147483648) %% 2147483646 + 1
  )

  ord <- order(registry$seed_token_sha256)
  used <- new.env(hash = TRUE, parent = emptyenv())
  final <- integer(nrow(registry))
  probes <- integer(nrow(registry))
  for (ii in ord) {
    s <- registry$raw_seed[ii]
    p <- 0L
    while (exists(as.character(s), envir = used, inherits = FALSE)) {
      s <- if (s == 2147483646L) 1L else s + 1L
      p <- p + 1L
    }
    assign(as.character(s), TRUE, envir = used)
    final[ii] <- s
    probes[ii] <- p
  }
  registry$seed <- final
  registry$collision_probes <- probes
  stopifnot(
    !anyDuplicated(registry$seed),
    !anyDuplicated(registry$seed_token_sha256),
    all(registry$seed > 0L)
  )
  registry
}

s5_build_seed_registry <- function(structural, rerun_cfg) {
  namespace <- as.character(rerun_cfg$seed_namespace)
  glabel <- as.character(rerun_cfg$global_seed_label)
  reps <- as.integer(rerun_cfg$replications_per_structural_cell)
  max_attempts <- as.integer(rerun_cfg$rng$max_bank_attempts)
  z <- vector("list", nrow(structural) * (2L + reps * (max_attempts + 2L)) + 1L)
  q <- 0L
  add <- function(sid, rep, stream, attempt = NA_integer_) {
    q <<- q + 1L
    z[[q]] <<- data.frame(
      structural_cell_id_full = sid, replication = rep,
      stream = stream, attempt = attempt,
      seed_token = s5_seed_token(namespace, glabel, sid, rep, stream, attempt),
      stringsAsFactors = FALSE
    )
  }
  for (h in seq_len(nrow(structural))) {
    sid <- structural$structural_cell_id_full[h]
    add(sid, NA_integer_, "quadrature_calibration")
    add(sid, NA_integer_, "quadrature_evaluation")
    for (r in seq_len(reps)) {
      for (a in seq_len(max_attempts)) add(sid, r, "item_bank", a)
      add(sid, r, "latent_abilities")
      add(sid, r, "response_uniforms")
    }
  }

  add("ALL_STRUCTURAL_CELLS", NA_integer_, "bootstrap_joint_form_blocks")
  out <- do.call(rbind, z[seq_len(q)])
  s5_allocate_content_hash_seeds(out)
}

s5_seed_lookup <- function(seed_registry, structural_cell_id_full,
                           replication = NA_integer_, stream,
                           attempt = NA_integer_) {
  use <- seed_registry$structural_cell_id_full == structural_cell_id_full &
    seed_registry$stream == stream &
    (is.na(seed_registry$replication) & is.na(replication) |
      seed_registry$replication == replication) &
    (is.na(seed_registry$attempt) & is.na(attempt) |
      seed_registry$attempt == attempt)
  use[is.na(use)] <- FALSE
  if (sum(use) != 1L) {
    stop("seed lookup did not identify one stream: ", stream,
      ", replication=", replication, ", attempt=", attempt,
      call. = FALSE
    )
  }
  seed_registry[use, , drop = FALSE]
}

s5_draw_item_bank <- function(args, seed_registry, structural_cell_id_full,
                              replication, max_tries = 20L) {
  failures <- character(0)
  for (a in seq_len(max_tries)) {
    sr <- s5_seed_lookup(
      seed_registry, structural_cell_id_full,
      replication, "item_bank", a
    )
    out <- tryCatch(
      do.call(
        IRTsimrel::sim_item_params,
        c(args, list(seed = as.integer(sr$seed)))
      ),
      error = function(e) e
    )
    if (!inherits(out, "error")) {
      return(list(
        value = out, status = "ok", attempts = a,
        seed = sr$seed, seed_token_sha256 = sr$seed_token_sha256,
        prior_failure_messages = paste(failures, collapse = " || ")
      ))
    }
    msg <- conditionMessage(out)
    failures <- c(failures, msg)
    if (!grepl(POOL_NA_MESSAGE, msg, fixed = TRUE)) stop(out)
  }
  stop("item bank failed all ", max_tries, " named attempts for form ",
    structural_cell_id_full, "/", replication,
    call. = FALSE
  )
}

s5_safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) stats::median(x) else NA_real_
}

s5_failure_codes <- function(status, p) {
  bad <- status != "ok" | !is.finite(p)
  bad[is.na(bad)] <- TRUE
  if (!any(bad)) {
    return("")
  }
  code <- as.character(status[bad])
  code[is.na(code) | !nzchar(code)] <- "missing_status"
  code[code == "ok" & !is.finite(p[bad])] <- "ok_nonfinite_p"
  tab <- table(code)
  paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = ";")
}

s5_arm_target <- function(base_cfg, arm) {
  x <- base_cfg$arms[[arm]]
  if (identical(x$kind, "uncontrolled")) NA_real_ else as.numeric(x$target_rho)
}

s5_aggregate_form <- function(res, row, condition_row, base_cfg, rerun_cfg,
                              replication, form_cluster_id_full,
                              variant_form_id_full, dif_items, cvals, realized,
                              bank_draw, ability_seed, uniform_seed, rng_kind) {
  arms <- as.character(unlist(rerun_cfg$arms))
  n_items <- as.integer(row$n_items)
  expected <- c(
    null_clean = n_items, power = length(dif_items),
    dif_form_clean_items = n_items - length(dif_items)
  )
  select_kind <- function(arm, kind) {
    if (is.null(res) || !nrow(res)) {
      return(res[FALSE, , drop = FALSE])
    }
    z <- res[res$arm == arm, , drop = FALSE]
    if (kind == "null_clean") {
      z[z$condition == "null", , drop = FALSE]
    } else if (kind == "power") {
      z[z$condition == "dif" & z$is_dif_item, , drop = FALSE]
    } else {
      z[z$condition == "dif" & !z$is_dif_item, , drop = FALSE]
    }
  }
  out <- vector("list", length(arms) * length(expected))
  q <- 0L
  for (arm in arms) {
    for (kind in names(expected)) {
      q <- q + 1L
      g <- select_kind(arm, kind)
      effect_g <- if (!is.null(res) && nrow(res)) {
        res[res$arm == arm & res$condition == "dif" & res$is_dif_item, ,
          drop = FALSE
        ]
      } else {
        res
      }
      reached <- arm %in% names(cvals) && is.finite(cvals[[arm]])
      ok <- if (nrow(g)) g$mh_status == "ok" & is.finite(g$mh_p) else logical(0)
      ok[is.na(ok)] <- FALSE
      n_attempt <- nrow(g)
      n_success <- sum(ok)
      n_reject <- if (n_success) sum(g$mh_p[ok] < rerun_cfg$alpha) else 0L
      status <- if (!reached) {
        "arm_unreachable"
      } else if (!n_attempt) {
        "missing_item_rows"
      } else if (!n_success) "all_mh_failed" else if (n_success < n_attempt) "partial_mh_failure" else "ok"
      qrow <- data.frame(
        schema_version = as.integer(rerun_cfg$schema_version),
        structural_cell_id_full = row$structural_cell_id_full,
        condition_id_full = condition_row$condition_id_full,
        source_cond_id = condition_row$source_cond_id,
        replication = as.integer(replication),
        form_cluster_id_full = form_cluster_id_full,
        variant_form_id_full = variant_form_id_full,
        variant = condition_row$variant,
        model = row$model, pool = row$pool,
        latent_shape = row$latent_shape, n_items = n_items,
        n_per_group = as.integer(row$n_per_group),
        arm = arm,
        arm_kind = as.character(base_cfg$arms[[arm]]$kind),
        target_rho = s5_arm_target(base_cfg, arm),
        arm_status = if (reached) "ok" else "target_unreachable",
        item_kind = kind,
        condition = if (kind == "null_clean") "null" else "dif",
        n_items_expected = as.integer(expected[[kind]]),
        n_mh_attempted = as.integer(n_attempt),
        n_mh_successful = as.integer(n_success),
        n_mh_failed = as.integer(n_attempt - n_success),
        n_mh_rejected = as.integer(n_reject),
        rate_mh_successful_only = if (n_success) n_reject / n_success else NA_real_,
        rate_mh_failures_as_nonrejections =
          if (n_attempt) n_reject / n_attempt else NA_real_,
        mh_form_status = status,
        mh_failure_codes = if (nrow(g)) {
          s5_failure_codes(g$mh_status, g$mh_p)
        } else {
          ""
        },
        c_value = if (reached) as.numeric(cvals[[arm]]) else NA_real_,
        realized_rho = if (arm %in% names(realized)) {
          as.numeric(realized[[arm]])
        } else {
          NA_real_
        },
        dif_scale = condition_row$dif_scale,
        purify = condition_row$purify,
        impact_mean_focal = condition_row$impact_mean_focal,
        effect_definition_id = condition_row$effect_definition_id,
        effect_scale_contract = condition_row$effect_scale_contract,
        effect_interpretation = condition_row$effect_interpretation,
        planted_theta_shift_median = if (!is.null(effect_g) && nrow(effect_g)) {
          s5_safe_median(effect_g$planted_shift)
        } else {
          NA_real_
        },
        planted_linear_predictor_shift =
          if (!is.null(effect_g) && nrow(effect_g)) {
            s5_safe_median(effect_g$planted_linear_predictor_shift)
          } else {
            NA_real_
          },
        planted_conditional_logodds =
          if (!is.null(effect_g) && nrow(effect_g)) {
            s5_safe_median(effect_g$planted_conditional_logodds)
          } else {
            NA_real_
          },
        planted_ets_delta_2pl_only =
          if (!is.null(effect_g) && nrow(effect_g)) {
            s5_safe_median(effect_g$planted_ets_delta)
          } else {
            NA_real_
          },
        n_anchor = if (nrow(g)) s5_safe_median(g$n_anchor) else NA_real_,
        bank_status = bank_draw$status,
        bank_attempts = as.integer(bank_draw$attempts),
        bank_retries = as.integer(bank_draw$attempts - 1L),
        bank_seed = as.integer(bank_draw$seed),
        bank_seed_token_sha256 = bank_draw$seed_token_sha256,
        ability_seed = as.integer(ability_seed$seed),
        ability_seed_token_sha256 = ability_seed$seed_token_sha256,
        uniform_seed = as.integer(uniform_seed$seed),
        uniform_seed_token_sha256 = uniform_seed$seed_token_sha256,
        rng_kind = rng_kind[1], rng_normal_kind = rng_kind[2],
        rng_sample_kind = rng_kind[3],
        uncertainty_scope = "MH_primary_form_cluster_rerun",
        logistic_scope = "retained_descriptive_points_only_no_form_cluster_inference",
        stringsAsFactors = FALSE
      )
      out[[q]] <- qrow
    }
  }
  do.call(rbind, out)
}

s5_validate_form_records <- function(x, catalog, rerun_cfg) {
  reps <- as.integer(rerun_cfg$replications_per_structural_cell)
  arms <- as.character(unlist(rerun_cfg$arms))
  variants <- names(rerun_cfg$variants)
  kinds <- c("null_clean", "power", "dif_form_clean_items")
  n_struct <- length(unique(catalog$structural_cell_id_full))
  expected_rows <- n_struct * reps * length(variants) * length(arms) *
    length(kinds)
  if (nrow(x) != expected_rows) {
    stop("form-record row count: expected ",
      expected_rows, ", observed ", nrow(x),
      call. = FALSE
    )
  }
  key <- interaction(x$condition_id_full, x$replication, x$arm, x$item_kind,
    drop = TRUE
  )
  if (anyDuplicated(key)) {
    stop("duplicate condition/form/arm/kind record",
      call. = FALSE
    )
  }
  if (!setequal(unique(x$condition_id_full), catalog$condition_id_full)) {
    stop("full condition-ID grid is incomplete", call. = FALSE)
  }
  if (!setequal(unique(x$variant), variants) ||
    !setequal(unique(x$arm), arms) || !setequal(unique(x$item_kind), kinds)) {
    stop("variant/arm/item-kind grid is incomplete", call. = FALSE)
  }
  obs <- table(x$structural_cell_id_full, x$variant, x$arm, x$item_kind)
  if (any(obs != reps)) {
    stop("planned/observed forms per cell are not ", reps,
      call. = FALSE
    )
  }
  reached <- x$arm_status == "ok"
  if (any(x$n_mh_attempted[reached] != x$n_items_expected[reached])) {
    stop("an reached form does not retain every attempted MH item", call. = FALSE)
  }
  if (any(x$n_mh_successful + x$n_mh_failed != x$n_mh_attempted)) {
    stop("MH successes and failures do not add to attempts", call. = FALSE)
  }
  if (any(x$n_mh_rejected > x$n_mh_successful)) {
    stop("MH rejections exceed successful tests", call. = FALSE)
  }
  if (!identical(sort(unique(x$impact_mean_focal[x$variant == "impact"])), 0.5)) {
    stop("impact sign/value gate failed", call. = FALSE)
  }
  if (any(x$impact_mean_focal[x$variant != "impact"] != 0)) {
    stop("non-impact variant carries impact", call. = FALSE)
  }
  lp <- x[x$variant == "logodds" & x$item_kind == "power", ]
  if (any(abs(lp$planted_linear_predictor_shift - 0.5) > 1e-10,
    na.rm = TRUE
  )) {
    stop("fixed-LP effect gate failed", call. = FALSE)
  }
  if (any(is.finite(lp$planted_conditional_logodds[lp$model != "2pl"]))) {
    stop("3PL rows claim a fixed conditional log-odds effect", call. = FALSE)
  }
  invisible(TRUE)
}

s5_rate_denominator <- function(x, estimand) {
  if (estimand == "successful_only") {
    x$n_mh_successful
  } else if (estimand == "failures_as_nonrejections") {
    x$n_mh_attempted
  } else {
    stop("unknown estimand: ", estimand, call. = FALSE)
  }
}

s5_stratified_ratio_se <- function(y, n, stratum) {
  keep <- is.finite(y) & is.finite(n) & n >= 0
  y <- y[keep]
  n <- n[keep]
  stratum <- stratum[keep]
  p <- sum(y) / sum(n)
  u <- y - p * n
  vv <- split(u, stratum)
  vtot <- sum(vapply(vv, function(z) {
    if (length(z) > 1L) length(z) * stats::var(z) else 0
  }, numeric(1)))
  c(estimate = p, mcse = sqrt(vtot) / sum(n))
}

s5_stratified_mean_se <- function(value, stratum) {
  keep <- is.finite(value)
  value <- value[keep]
  stratum <- stratum[keep]
  z <- split(value, stratum)
  n <- length(value)
  vtot <- sum(vapply(z, function(q) {
    if (length(q) > 1L) length(q) * stats::var(q) else 0
  }, numeric(1)))
  c(estimate = mean(value), mcse = sqrt(vtot) / n)
}

s5_cluster_rate_table <- function(x, ci_level = 0.95) {
  zcrit <- stats::qnorm(1 - (1 - ci_level) / 2)
  groups <- split(x, list(x$variant, x$arm, x$item_kind), drop = TRUE)
  out <- list()
  q <- 0L
  for (g in groups) {
    for (estimand in c(
      "successful_only",
      "failures_as_nonrejections"
    )) {
      q <- q + 1L
      den <- s5_rate_denominator(g, estimand)
      ss <- s5_stratified_ratio_se(
        g$n_mh_rejected, den,
        g$structural_cell_id_full
      )
      out[[q]] <- data.frame(
        variant = g$variant[1], arm = g$arm[1], item_kind = g$item_kind[1],
        estimand = estimand,
        n_structural_cells = length(unique(g$structural_cell_id_full)),
        n_form_clusters_planned = length(unique(g$form_cluster_id_full)),
        n_form_clusters_contributing = sum(den > 0),
        n_form_clusters_unreachable = sum(g$arm_status != "ok"),
        n_form_clusters_with_mh_failure = sum(g$n_mh_failed > 0),
        n_mh_rejected = sum(g$n_mh_rejected),
        n_mh_successful = sum(g$n_mh_successful),
        n_mh_failed = sum(g$n_mh_failed),
        n_mh_attempted = sum(g$n_mh_attempted),
        mh_failure_rate = sum(g$n_mh_failed) / sum(g$n_mh_attempted),
        estimate = unname(ss["estimate"]),
        mcse_form_cluster_sandwich = unname(ss["mcse"]),
        ci_low_form_cluster_normal = max(0, ss["estimate"] - zcrit * ss["mcse"]),
        ci_high_form_cluster_normal = min(1, ss["estimate"] + zcrit * ss["mcse"]),
        realized_rho_median = stats::median(g$realized_rho, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, out)
}

s5_paired_ratio_influence <- function(a, b, estimand) {
  na <- s5_rate_denominator(a, estimand)
  nb <- s5_rate_denominator(b, estimand)
  pa <- sum(a$n_mh_rejected) / sum(na)
  pb <- sum(b$n_mh_rejected) / sum(nb)
  infl <- (b$n_mh_rejected - pb * nb) / sum(nb) -
    (a$n_mh_rejected - pa * na) / sum(na)
  list(estimate = pb - pa, influence = infl)
}

s5_arm_covariance_table <- function(x) {
  out <- list()
  q <- 0L
  for (v in unique(x$variant)) {
    for (k in unique(x$item_kind)) {
      for (estimand in c("successful_only", "failures_as_nonrejections")) {
        g <- x[x$variant == v & x$item_kind == k, ]
        arms <- sort(unique(g$arm))
        for (ij in utils::combn(arms, 2L, simplify = FALSE)) {
          a <- g[g$arm == ij[1], ]
          b <- g[g$arm == ij[2], ]
          m <- merge(a, b, by = c(
            "structural_cell_id_full", "replication",
            "form_cluster_id_full"
          ), suffixes = c("_a", "_b"))
          da <- if (estimand == "successful_only") {
            m$n_mh_successful_a
          } else {
            m$n_mh_attempted_a
          }
          db <- if (estimand == "successful_only") {
            m$n_mh_successful_b
          } else {
            m$n_mh_attempted_b
          }
          ra <- m$n_mh_rejected_a / da
          rb <- m$n_mh_rejected_b / db
          keep <- is.finite(ra) & is.finite(rb)
          strata <- split(seq_len(nrow(m))[keep], m$structural_cell_id_full[keep])
          covs <- vapply(strata, function(ii) {
            if (length(ii) > 1L) {
              stats::cov(ra[ii], rb[ii])
            } else {
              0
            }
          }, numeric(1))
          vars_a <- vapply(strata, function(ii) {
            if (length(ii) > 1L) {
              stats::var(ra[ii])
            } else {
              0
            }
          }, numeric(1))
          vars_b <- vapply(strata, function(ii) {
            if (length(ii) > 1L) {
              stats::var(rb[ii])
            } else {
              0
            }
          }, numeric(1))
          nh <- vapply(strata, length, integer(1))
          N <- sum(nh)
          q <- q + 1L
          out[[q]] <- data.frame(
            variant = v, item_kind = k, estimand = estimand,
            arm_a = ij[1], arm_b = ij[2], n_paired_forms = N,
            mean_within_cell_covariance = weighted.mean(covs, nh),
            covariance_of_arm_rate_estimators = sum(nh * covs) / N^2,
            correlation_of_arm_rate_estimators = sum(nh * covs) /
              sqrt(sum(nh * vars_a) * sum(nh * vars_b)),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  do.call(rbind, out)
}

s5_paired_arm_contrast_table <- function(x, reference = "uncontrolled",
                                         ci_level = 0.95) {
  zcrit <- stats::qnorm(1 - (1 - ci_level) / 2)
  out <- list()
  q <- 0L
  for (v in unique(x$variant)) {
    for (k in unique(x$item_kind)) {
      g <- x[x$variant == v & x$item_kind == k, ]
      for (arm in setdiff(unique(g$arm), reference)) {
        for (estimand in c("successful_only", "failures_as_nonrejections")) {
          a <- g[g$arm == reference, ]
          b <- g[g$arm == arm, ]
          m <- merge(a, b, by = c(
            "structural_cell_id_full", "replication",
            "form_cluster_id_full"
          ), suffixes = c("_a", "_b"))
          na <- if (estimand == "successful_only") {
            m$n_mh_successful_a
          } else {
            m$n_mh_attempted_a
          }
          nb <- if (estimand == "successful_only") {
            m$n_mh_successful_b
          } else {
            m$n_mh_attempted_b
          }
          pa <- sum(m$n_mh_rejected_a) / sum(na)
          pb <- sum(m$n_mh_rejected_b) / sum(nb)
          infl <- (m$n_mh_rejected_b - pb * nb) / sum(nb) -
            (m$n_mh_rejected_a - pa * na) / sum(na)
          strata <- split(infl, m$structural_cell_id_full)
          vv <- sum(vapply(strata, function(z) {
            if (length(z) > 1L) {
              length(z) * stats::var(z)
            } else {
              0
            }
          }, numeric(1)))
          se <- sqrt(vv)
          delta <- pb - pa
          q <- q + 1L
          out[[q]] <- data.frame(
            variant = v, item_kind = k, estimand = estimand,
            reference_arm = reference, comparison_arm = arm,
            n_joint_form_records = nrow(m),
            n_reference_contributing_forms = sum(na > 0),
            n_comparison_contributing_forms = sum(nb > 0),
            reference_rate = pa,
            comparison_rate = pb, paired_difference = delta,
            mcse_paired_form_cluster = se,
            ci_low_paired_normal = delta - zcrit * se,
            ci_high_paired_normal = delta + zcrit * se,
            pairing = "same_structural_cell_form_bank_abilities_uniforms",
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  do.call(rbind, out)
}

s5_bootstrap_joint <- function(x, structural, rerun_cfg, seed_registry) {
  B <- as.integer(rerun_cfg$inference$bootstrap_replications)
  variants <- names(rerun_cfg$variants)
  arms <- as.character(unlist(rerun_cfg$arms))
  kinds <- c("null_clean", "power", "dif_form_clean_items")
  estimands <- c("successful_only", "failures_as_nonrejections")
  reps <- as.integer(rerun_cfg$replications_per_structural_cell)
  H <- nrow(structural)
  Q <- H * reps
  hmap <- setNames(seq_len(H), structural$structural_cell_id_full)
  vmap <- setNames(seq_along(variants), variants)
  amap <- setNames(seq_along(arms), arms)
  kmap <- setNames(seq_along(kinds), kinds)
  M <- length(variants) * length(arms) * length(kinds)
  Y <- matrix(NA_real_, Q, M)
  DS <- Y
  DA <- Y
  col_index <- function(v, a, k) {
    (vmap[[v]] - 1L) * length(arms) * length(kinds) +
      (amap[[a]] - 1L) * length(kinds) + kmap[[k]]
  }
  for (i in seq_len(nrow(x))) {
    h <- hmap[[x$structural_cell_id_full[i]]]
    qi <- (h - 1L) * reps + x$replication[i]
    mi <- col_index(x$variant[i], x$arm[i], x$item_kind[i])
    Y[qi, mi] <- x$n_mh_rejected[i]
    DS[qi, mi] <- x$n_mh_successful[i]
    DA[qi, mi] <- x$n_mh_attempted[i]
  }
  stopifnot(!anyNA(Y), !anyNA(DS), !anyNA(DA))

  meta <- structural[rep(seq_len(H), each = reps), , drop = FALSE]
  rate_template <- do.call(rbind, lapply(estimands, function(e) {
    do.call(rbind, lapply(variants, function(v) {
      do.call(rbind, lapply(arms, function(a) {
        data.frame(
          variant = v, arm = a, item_kind = kinds,
          estimand = e, stringsAsFactors = FALSE
        )
      }))
    }))
  }))
  rate_draws <- rate_template[rep(seq_len(nrow(rate_template)), times = B), ]
  rate_draws$bootstrap_replication <- rep(seq_len(B), each = nrow(rate_template))
  rate_draws$estimate <- NA_real_
  factor_names <- c(
    length = "n_items", model = "model",
    sample_size = "n_per_group", pool = "pool",
    shape = "latent_shape"
  )
  high <- c(
    length = "40", model = "3pl_g20", sample_size = "1000",
    pool = "irw", shape = "skew_pos"
  )
  low <- c(
    length = "20", model = "2pl", sample_size = "500",
    pool = "parametric", shape = "normal"
  )
  factor_template <- do.call(rbind, lapply(estimands, function(e) {
    do.call(rbind, lapply(variants, function(v) {
      data.frame(
        variant = v, contrast = names(factor_names),
        factor = unname(factor_names), level_high = unname(high),
        level_low = unname(low), estimand = e,
        stringsAsFactors = FALSE
      )
    }))
  }))
  factor_draws <- factor_template[
    rep(seq_len(nrow(factor_template)), times = B),
  ]
  factor_draws$bootstrap_replication <-
    rep(seq_len(B), each = nrow(factor_template))
  factor_draws$delta_uncontrolled <- NA_real_
  factor_draws$delta_targeted <- NA_real_
  factor_draws$did <- NA_real_
  sr <- s5_seed_lookup(
    seed_registry, "ALL_STRUCTURAL_CELLS", NA_integer_,
    "bootstrap_joint_form_blocks", NA_integer_
  )
  RNGkind(
    as.character(rerun_cfg$rng$kind),
    as.character(rerun_cfg$rng$normal_kind),
    as.character(rerun_cfg$rng$sample_kind)
  )
  set.seed(as.integer(sr$seed))
  for (b in seq_len(B)) {
    rate_offset <- (b - 1L) * nrow(rate_template)
    factor_offset <- (b - 1L) * nrow(factor_template)
    w <- unlist(
      lapply(seq_len(H), function(h) {
        tabulate(sample.int(reps, reps, replace = TRUE), nbins = reps)
      }),
      use.names = FALSE
    )
    wy <- as.numeric(crossprod(w, Y))
    for (e in seq_along(estimands)) {
      denmat <- if (e == 1L) DS else DA
      wd <- as.numeric(crossprod(w, denmat))
      rates <- wy / wd
      ridx <- rate_offset + (e - 1L) * M + seq_len(M)
      rate_draws$estimate[ridx] <- rates

      ref_arms <- c("uncontrolled", as.character(rerun_cfg$reference_arm))
      for (vi in seq_along(variants)) {
        for (fi in seq_along(factor_names)) {
          v <- variants[vi]
          fn <- names(factor_names)[fi]
          fac <- factor_names[[fn]]
          vals <- list()
          for (a in ref_arms) {
            mi <- col_index(v, a, "power")
            den <- denmat[, mi]
            rform <- Y[, mi] / den
            vals[[a]] <- rform
          }
          pair_ok <- is.finite(vals[[ref_arms[1]]]) &
            is.finite(vals[[ref_arms[2]]])
          lev <- as.character(meta[[fac]])
          calc <- function(z, level) {
            use <- lev == level & pair_ok
            sum(w[use] * z[use]) / sum(w[use])
          }
          du <- calc(vals[[ref_arms[1]]], high[[fn]]) -
            calc(vals[[ref_arms[1]]], low[[fn]])
          dt <- calc(vals[[ref_arms[2]]], high[[fn]]) -
            calc(vals[[ref_arms[2]]], low[[fn]])
          ti <- (e - 1L) * length(variants) * length(factor_names) +
            (vi - 1L) * length(factor_names) + fi
          fidx <- factor_offset + ti
          factor_draws$delta_uncontrolled[fidx] <- du
          factor_draws$delta_targeted[fidx] <- dt
          factor_draws$did[fidx] <- dt - du
        }
      }
    }
  }
  list(
    rate_draws = rate_draws,
    factor_draws = factor_draws,
    bootstrap_seed = sr$seed,
    bootstrap_seed_token_sha256 = sr$seed_token_sha256,
    bootstrap_replications = B,
    resampling_unit = "form_cluster_block_within_structural_cell",
    covariance_preserved = "all_arms_and_all_variants"
  )
}

s5_add_bootstrap_rate_intervals <- function(rate_table, rate_draws,
                                            ci_level = 0.95) {
  alpha <- (1 - ci_level) / 2
  bs <- do.call(rbind, lapply(split(rate_draws,
    list(
      rate_draws$variant, rate_draws$arm, rate_draws$item_kind,
      rate_draws$estimand
    ),
    drop = TRUE
  ), function(g) {
    data.frame(
      variant = g$variant[1], arm = g$arm[1], item_kind = g$item_kind[1],
      estimand = g$estimand[1],
      mcse_cell_stratified_bootstrap = stats::sd(g$estimate),
      ci_low_cell_stratified_bootstrap = stats::quantile(g$estimate, alpha),
      ci_high_cell_stratified_bootstrap = stats::quantile(g$estimate, 1 - alpha),
      stringsAsFactors = FALSE
    )
  }))
  merge(rate_table, bs,
    by = c("variant", "arm", "item_kind", "estimand"),
    sort = FALSE
  )
}

s5_factor_point_table <- function(x, reference_arm = "target_075") {
  power <- x[x$item_kind == "power" &
    x$arm %in% c("uncontrolled", reference_arm), ]
  out <- list()
  q <- 0L
  specs <- list(
    length = c(factor = "n_items", high = "40", low = "20"),
    model = c(factor = "model", high = "3pl_g20", low = "2pl"),
    sample_size = c(factor = "n_per_group", high = "1000", low = "500"),
    pool = c(factor = "pool", high = "irw", low = "parametric"),
    shape = c(factor = "latent_shape", high = "skew_pos", low = "normal")
  )
  for (v in unique(power$variant)) {
    for (nm in names(specs)) {
      for (estimand in c("successful_only", "failures_as_nonrejections")) {
        g <- power[power$variant == v, ]
        by <- c(
          "structural_cell_id_full", "replication", "form_cluster_id_full",
          "model", "pool", "latent_shape", "n_items", "n_per_group"
        )
        u <- g[g$arm == "uncontrolled", ]
        t <- g[g$arm == reference_arm, ]
        m <- merge(u, t, by = by, suffixes = c("_u", "_t"))
        den_u <- if (estimand == "successful_only") {
          m$n_mh_successful_u
        } else {
          m$n_mh_attempted_u
        }
        den_t <- if (estimand == "successful_only") {
          m$n_mh_successful_t
        } else {
          m$n_mh_attempted_t
        }
        m$r_u <- m$n_mh_rejected_u / den_u
        m$r_t <- m$n_mh_rejected_t / den_t
        m <- m[is.finite(m$r_u) & is.finite(m$r_t), , drop = FALSE]
        s <- specs[[nm]]
        lev <- as.character(m[[s[["factor"]]]])
        avg <- function(z, level) mean(z[lev == level])
        du <- avg(m$r_u, s[["high"]]) - avg(m$r_u, s[["low"]])
        dt <- avg(m$r_t, s[["high"]]) - avg(m$r_t, s[["low"]])
        q <- q + 1L
        out[[q]] <- data.frame(
          contrast = nm, variant = v, factor = s[["factor"]],
          level_high = s[["high"]], level_low = s[["low"]], estimand = estimand,
          n_complete_paired_forms = nrow(m),
          delta_uncontrolled = du, delta_targeted = dt, did = dt - du,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, out)
}

s5_add_factor_bootstrap <- function(point, draws, ci_level = 0.95) {
  alpha <- (1 - ci_level) / 2
  summ <- do.call(rbind, lapply(split(draws,
    list(draws$variant, draws$contrast, draws$estimand),
    drop = TRUE
  ), function(g) {
    one <- function(v) {
      c(
        mcse = stats::sd(v),
        low = unname(stats::quantile(v, alpha)),
        high = unname(stats::quantile(v, 1 - alpha))
      )
    }
    u <- one(g$delta_uncontrolled)
    t <- one(g$delta_targeted)
    d <- one(g$did)
    data.frame(
      variant = g$variant[1], contrast = g$contrast[1],
      estimand = g$estimand[1],
      delta_uncontrolled_mcse = u["mcse"], delta_uncontrolled_ci_low = u["low"],
      delta_uncontrolled_ci_high = u["high"],
      delta_targeted_mcse = t["mcse"], delta_targeted_ci_low = t["low"],
      delta_targeted_ci_high = t["high"],
      did_mcse = d["mcse"], did_ci_low = d["low"], did_ci_high = d["high"],
      stringsAsFactors = FALSE
    )
  }))
  merge(point, summ, by = c("variant", "contrast", "estimand"), sort = FALSE)
}

s5_variant_contrast_table <- function(factor_table, factor_draws,
                                      baseline = "baseline", ci_level = 0.95) {
  alpha <- (1 - ci_level) / 2
  alternatives <- setdiff(unique(factor_table$variant), baseline)
  out <- list()
  q <- 0L
  for (alt in alternatives) {
    for (co in unique(factor_table$contrast)) {
      for (est in unique(factor_table$estimand)) {
        a <- factor_table[factor_table$variant == alt &
          factor_table$contrast == co & factor_table$estimand == est, ]
        b <- factor_table[factor_table$variant == baseline &
          factor_table$contrast == co & factor_table$estimand == est, ]
        da <- factor_draws[factor_draws$variant == alt &
          factor_draws$contrast == co & factor_draws$estimand == est, ]
        db <- factor_draws[factor_draws$variant == baseline &
          factor_draws$contrast == co & factor_draws$estimand == est, ]
        m <- merge(da, db, by = "bootstrap_replication", suffixes = c("_alt", "_base"))
        vd <- m$did_alt - m$did_base
        q <- q + 1L
        out[[q]] <- data.frame(
          contrast = co, alternative_variant = alt, baseline_variant = baseline,
          estimand = est, did_alternative = a$did, did_baseline = b$did,
          paired_variant_did_difference = a$did - b$did,
          mcse_joint_cell_stratified_bootstrap = stats::sd(vd),
          ci_low_joint_cell_stratified_bootstrap = stats::quantile(vd, alpha),
          ci_high_joint_cell_stratified_bootstrap = stats::quantile(vd, 1 - alpha),
          pairing = "same_form_bank_abilities_uniforms_and_bootstrap_blocks",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, out)
}
