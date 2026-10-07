# Paired response generation and Mantel-Haenszel and logistic DIF tests.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

REST_SCORE_BREAKS <- function(n_items) {
  b <- unique(c(-1, seq(2, n_items - 1, by = 3), n_items))
  if (length(b) < 3L) b <- c(-1, floor(n_items / 2), n_items)
  b
}

irt_prob <- function(theta, beta, lambda, guessing = NULL) {
  stopifnot(length(beta) == length(lambda))
  if (length(lambda) == 1L) lambda <- rep(lambda, length(beta))
  g <- if (is.null(guessing)) rep(0, length(beta)) else if (length(guessing) == 1L) rep(guessing, length(beta)) else guessing
  stopifnot(length(g) == length(beta), all(g >= 0), all(g < 1))
  eta <- outer(theta, beta, "-") * rep(lambda, each = length(theta))
  rep(g, each = length(theta)) +
    rep(1 - g, each = length(theta)) * stats::plogis(eta)
}

mh_rest_score <- function(resp_ref, resp_foc, item, breaks = NULL,
                          anchor = NULL) {
  n_items <- ncol(resp_ref)
  use_cols <- if (is.null(anchor)) seq_len(n_items) else anchor
  use_cols <- setdiff(use_cols, item)
  if (!length(use_cols)) {
    return(list(
      status = "empty_matching_set", statistic = NA_real_,
      p_value = NA_real_, hand_statistic = NA_real_,
      hand_p_value = NA_real_, informative_strata = 0L,
      excluded_strata = NA_integer_
    ))
  }

  if (is.null(breaks)) breaks <- REST_SCORE_BREAKS(length(use_cols) + 1L)
  labs <- paste0("s", seq_len(length(breaks) - 1L))

  rest <- c(
    rowSums(resp_ref[, use_cols, drop = FALSE]),
    rowSums(resp_foc[, use_cols, drop = FALSE])
  )
  y <- c(resp_ref[, item], resp_foc[, item])
  grp <- factor(c(rep("reference", nrow(resp_ref)), rep("focal", nrow(resp_foc))),
    levels = c("reference", "focal")
  )
  strata <- cut(rest,
    breaks = breaks, labels = labs, include.lowest = TRUE,
    right = TRUE
  )

  arr <- array(0L,
    dim = c(2L, 2L, length(labs)),
    dimnames = list(
      group = c("reference", "focal"),
      response = c("0", "1"), stratum = labs
    )
  )
  for (k in seq_along(labs)) {
    use <- !is.na(strata) & strata == labs[k]
    if (any(use)) {
      arr[, , k] <- table(
        factor(grp[use], levels = c("reference", "focal")),
        factor(y[use], levels = c(0, 1))
      )
    }
  }
  informative <- vapply(seq_along(labs), function(k) {
    tab <- arr[, , k]
    sum(tab) > 1L && all(rowSums(tab) > 0L) && all(colSums(tab) > 0L)
  }, logical(1))

  if (!any(informative)) {
    return(list(
      status = "no_informative_strata", statistic = NA_real_,
      p_value = NA_real_, hand_statistic = NA_real_,
      hand_p_value = NA_real_, informative_strata = 0L,
      excluded_strata = length(labs)
    ))
  }

  inf_arr <- arr[, , informative, drop = FALSE]
  res <- tryCatch(stats::mantelhaen.test(inf_arr, correct = TRUE),
    error = function(e) NULL, warning = function(w) NULL
  )

  delta <- 0
  variance <- 0
  for (k in seq_len(dim(inf_arr)[3L])) {
    tab <- inf_arr[, , k]
    n <- sum(tab)
    delta <- delta + tab[1L, 1L] - rowSums(tab)[1L] * colSums(tab)[1L] / n
    variance <- variance + prod(rowSums(tab)) * prod(colSums(tab)) / (n^2 * (n - 1))
  }
  yates <- if (abs(delta) >= 0.5) 0.5 else 0
  hand_stat <- (abs(delta) - yates)^2 / variance
  hand_p <- stats::pchisq(hand_stat, df = 1, lower.tail = FALSE)

  list(
    status = if (is.null(res)) "test_failed" else "ok",
    statistic = if (is.null(res)) NA_real_ else unname(res$statistic),
    p_value = if (is.null(res)) NA_real_ else unname(res$p.value),
    hand_statistic = unname(hand_stat), hand_p_value = unname(hand_p),
    informative_strata = sum(informative),
    excluded_strata = sum(!informative)
  )
}

lr_dif <- function(resp_ref, resp_foc, item, matching = c("rest", "total")) {
  matching <- match.arg(matching)
  drop_item <- if (matching == "rest") -item else seq_len(ncol(resp_ref))
  score <- c(
    rowSums(resp_ref[, drop_item, drop = FALSE]),
    rowSums(resp_foc[, drop_item, drop = FALSE])
  )
  y <- c(resp_ref[, item], resp_foc[, item])
  g <- c(rep(0L, nrow(resp_ref)), rep(1L, nrow(resp_foc)))

  if (length(unique(y)) < 2L) {
    return(list(
      status = "constant_response", p_uniform = NA_real_,
      p_nonuniform = NA_real_, beta_group = NA_real_,
      beta_interaction = NA_real_
    ))
  }
  fit <- function(fml) {
    tryCatch(suppressWarnings(
      stats::glm(fml, family = stats::binomial())
    ), error = function(e) NULL)
  }
  m0 <- fit(y ~ score)
  m1 <- fit(y ~ score + g)
  m2 <- fit(y ~ score + g + score:g)
  if (is.null(m0) || is.null(m1) || is.null(m2) ||
    !m0$converged || !m1$converged || !m2$converged) {
    return(list(
      status = "fit_failed", p_uniform = NA_real_,
      p_nonuniform = NA_real_, beta_group = NA_real_,
      beta_interaction = NA_real_
    ))
  }
  lrt <- function(small, big, df = 1L) {
    stats::pchisq(small$deviance - big$deviance, df = df, lower.tail = FALSE)
  }

  cf1 <- stats::coef(m1)
  cf2 <- stats::coef(m2)
  list(
    status = "ok",
    p_uniform = lrt(m0, m1),
    p_nonuniform = lrt(m1, m2),
    beta_group = unname(cf1[["g"]]),
    beta_interaction = if ("score:g" %in% names(cf2)) {
      unname(cf2[["score:g"]])
    } else {
      NA_real_
    }
  )
}

dif_replicate <- function(bank, c_by_arm, dif_items, dif_shift = 0.5,
                          n_per_group = 500L, seed = 1L,
                          uniform_seed = NULL,
                          tests = c("mh", "lr"),
                          dif_scale = c("theta", "logodds"),
                          impact = 0, purify = FALSE, alpha = 0.05) {
  dif_scale <- match.arg(dif_scale)
  n_items <- length(bank$beta)
  set.seed(seed)
  theta_ref <- stats::rnorm(n_per_group)
  theta_foc <- stats::rnorm(n_per_group, mean = impact)
  if (is.null(uniform_seed)) uniform_seed <- seed + 1L
  set.seed(as.integer(uniform_seed))
  u_ref <- matrix(stats::runif(n_per_group * n_items), n_per_group)
  u_foc <- matrix(stats::runif(n_per_group * n_items), n_per_group)

  out <- list()
  for (arm in names(c_by_arm)) {
    cv <- c_by_arm[[arm]]
    lam <- bank$lambda_base * cv

    shift_vec <- if (dif_scale == "theta") {
      rep(dif_shift, n_items)
    } else {
      dif_shift / lam
    }
    planted_linear_predictor_shift <- mean((lam * shift_vec)[dif_items])
    is_2pl <- is.null(bank$guessing) || all(bank$guessing == 0)
    planted_conditional_logodds <- if (is_2pl) {
      planted_linear_predictor_shift
    } else {
      NA_real_
    }
    planted_ets_delta <- if (is_2pl) {
      2.35 * planted_conditional_logodds
    } else {
      NA_real_
    }
    p_ref <- irt_prob(theta_ref, bank$beta, lam, bank$guessing)
    r_ref <- 1L * (u_ref < p_ref)

    for (cond in c("null", "dif")) {
      beta_foc <- bank$beta
      if (cond == "dif") {
        beta_foc[dif_items] <- beta_foc[dif_items] + shift_vec[dif_items]
      }
      p_foc <- irt_prob(theta_foc, beta_foc, lam, bank$guessing)
      r_foc <- 1L * (u_foc < p_foc)

      anchor <- NULL
      if (purify && "mh" %in% tests) {
        first <- vapply(seq_len(n_items), function(i) {
          m <- mh_rest_score(r_ref, r_foc, i)
          if (is.na(m$p_value)) NA_real_ else m$p_value
        }, numeric(1))
        flagged <- which(!is.na(first) & first < alpha)
        anchor <- setdiff(seq_len(n_items), flagged)

        if (length(anchor) < 3L) anchor <- NULL
      }

      for (i in seq_len(n_items)) {
        row <- data.frame(
          arm = arm, condition = cond, item = i,
          is_dif_item = i %in% dif_items,
          c_value = cv, dif_scale = dif_scale,
          impact = impact, purified = purify,
          planted_shift = shift_vec[i],
          planted_linear_predictor_shift =
            planted_linear_predictor_shift,
          planted_conditional_logodds =
            planted_conditional_logodds,
          planted_ets_delta = planted_ets_delta,
          planted_effect_scale = if (is_2pl) {
            "conditional_log_odds_and_linear_predictor"
          } else {
            "linear_predictor_only"
          },
          n_anchor = if (is.null(anchor)) {
            NA_integer_
          } else {
            length(anchor)
          },
          stringsAsFactors = FALSE
        )
        if ("mh" %in% tests) {
          m <- mh_rest_score(r_ref, r_foc, i, anchor = anchor)
          row$mh_status <- m$status
          row$mh_p <- m$p_value
          row$mh_hand_p <- m$hand_p_value
          row$mh_strata <- m$informative_strata
        }
        if ("lr" %in% tests) {
          l <- lr_dif(r_ref, r_foc, i)
          row$lr_status <- l$status
          row$lr_p_uniform <- l$p_uniform
          row$lr_p_nonuniform <- l$p_nonuniform
        }
        out[[length(out) + 1L]] <- row
      }
    }
  }
  do.call(rbind, out)
}
