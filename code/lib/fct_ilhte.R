# Response generation, mixed models, and diagnostics for item-level treatment effects.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

# Draw one response matrix under item-specific treatment effects.
# Treatment is assigned to persons; each item has one effect shared by those persons.
ilhte_data <- function(n_sub = 500L, n_items = 20L, te_mean = 0, te_sd = 0,
                       disc = 1.7, seed = 1L) {
  set.seed(seed)
  theta <- stats::rnorm(n_sub)
  b <- stats::rnorm(n_items)
  tau <- stats::rnorm(n_items, te_mean, te_sd)
  treat <- rep(0L, n_sub)
  treat[sample.int(n_sub, floor(n_sub / 2))] <- 1L

  lin <- outer(theta, b, "-") * disc +
    outer(treat, tau, "*")
  y <- matrix(stats::rbinom(length(lin), 1L, stats::plogis(lin)), nrow = n_sub)

  data.frame(
    person = rep(seq_len(n_sub), times = n_items),
    item = rep(seq_len(n_items), each = n_sub),
    treat = rep(treat, times = n_items),
    correct = as.vector(y),
    stringsAsFactors = FALSE
  )
}

# Fit constant and varying item-effect models to the same response data.
ilhte_fit <- function(d) {
  ctrl <- lme4::glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
  fit_or_null <- function(fml) {
    tryCatch(
      suppressMessages(suppressWarnings(
        lme4::glmer(fml, data = d, family = stats::binomial, control = ctrl)
      )),
      error = function(e) NULL
    )
  }

  m_const <- fit_or_null(correct ~ treat + (1 | person) + (1 | item))
  m_vary <- fit_or_null(correct ~ treat + (1 | person) + (treat | item))

  grab <- function(m, label) {
    if (is.null(m)) {
      return(data.frame(
        method = label, estimate = NA_real_, std.error = NA_real_,
        statistic = NA_real_, p.value = NA_real_,
        converged = FALSE, optimizer_converged = FALSE,
        singular = NA, boundary_cholesky = NA,
        stringsAsFactors = FALSE
      ))
    }
    cf <- summary(m)$coefficients
    if (!"treat" %in% rownames(cf)) {
      return(data.frame(
        method = label, estimate = NA_real_, std.error = NA_real_,
        statistic = NA_real_, p.value = NA_real_,
        converged = FALSE, optimizer_converged = FALSE,
        singular = NA, boundary_cholesky = NA,
        stringsAsFactors = FALSE
      ))
    }
    opt_code <- m@optinfo$conv$opt
    optimizer_converged <- length(opt_code) > 0L &&
      all(is.finite(opt_code)) && all(opt_code == 0)
    lme4_messages <- as.character(m@optinfo$conv$lme4$messages)
    nonboundary_messages <- lme4_messages[
      !grepl("boundary \\(singular\\)", lme4_messages, ignore.case = TRUE)
    ]
    boundary_cholesky <- any(is.finite(lme4::getME(m, "lower")) &
      abs(lme4::getME(m, "theta") -
        lme4::getME(m, "lower")) <= 1e-7)
    data.frame(
      method = label, estimate = unname(cf["treat", 1]),
      std.error = unname(cf["treat", 2]),
      statistic = unname(cf["treat", 3]),
      p.value = unname(cf["treat", 4]),
      converged = optimizer_converged && !length(nonboundary_messages),
      optimizer_converged = optimizer_converged,
      singular = lme4::isSingular(m, tol = 1e-4),
      boundary_cholesky = boundary_cholesky,
      stringsAsFactors = FALSE
    )
  }

  lrt_stat <- if (is.null(m_const) || is.null(m_vary)) {
    NA_real_
  } else {
    tryCatch(
      max(0, 2 * (as.numeric(stats::logLik(m_vary)) -
        as.numeric(stats::logLik(m_const)))),
      error = function(e) NA_real_
    )
  }
  sig_het_df2 <- if (is.finite(lrt_stat)) {
    stats::pchisq(lrt_stat, df = 2, lower.tail = FALSE)
  } else {
    NA_real_
  }
  sig_het_chibar <- ilhte_chibar_p(lrt_stat)

  out <- rbind(grab(m_const, "Constant"), grab(m_vary, "Varying"))

  out$sig_het <- sig_het_df2
  out$sig_het_lrt <- lrt_stat
  out$sig_het_df2 <- sig_het_df2
  out$sig_het_chibar <- sig_het_chibar
  out
}

ilhte_chibar_p <- function(lrt_stat) {
  out <- rep(NA_real_, length(lrt_stat))
  ok <- is.finite(lrt_stat) & lrt_stat >= 0
  out[ok] <- 0.5 * stats::pchisq(lrt_stat[ok], df = 1, lower.tail = FALSE) +
    0.5 * stats::pchisq(lrt_stat[ok], df = 2, lower.tail = FALSE)
  out
}

ilhte_objective_diagnostics <- function(
  fn, par, lower = rep(-Inf, length(par)), upper = rep(Inf, length(par)),
  delta = 1e-4, boundary_tol = 1e-7, gradient_tol = 0.002,
  hessian_eigen_tol = 1e-6, hessian_ratio_tol = 1e-6
) {
  stopifnot(
    is.function(fn), length(par) >= 1L,
    length(lower) == length(par), length(upper) == length(par),
    delta > 0, boundary_tol >= 0, gradient_tol > 0
  )
  n <- length(par)
  f0 <- suppressWarnings(fn(par))
  if (length(f0) != 1L || !is.finite(f0)) {
    return(data.frame(
      objective = NA_real_, n_parameters = n, n_free = NA_integer_,
      n_active_lower = NA_integer_, n_active_upper = NA_integer_,
      max_abs_gradient = NA_real_, max_abs_scaled_gradient = NA_real_,
      max_min_gradient = NA_real_, gradient_ok = FALSE,
      kkt_min_margin = NA_real_, kkt_ok = FALSE,
      hessian_min_eigen = NA_real_, hessian_max_eigen = NA_real_,
      hessian_eigen_ratio = NA_real_, hessian_ok = FALSE,
      numerical_ok = FALSE, stringsAsFactors = FALSE
    ))
  }

  active_lower <- is.finite(lower) & par <= lower + boundary_tol
  active_upper <- is.finite(upper) & par >= upper - boundary_tol
  free <- !(active_lower | active_upper)
  free_idx <- which(free)
  gradient <- rep(NA_real_, n)

  for (j in seq_len(n)) {
    step <- rep(0, n)
    step[j] <- delta
    gradient[j] <- if (active_lower[j]) {
      (-3 * f0 + 4 * fn(par + step) - fn(par + 2 * step)) / (2 * delta)
    } else if (active_upper[j]) {
      (3 * f0 - 4 * fn(par - step) + fn(par - 2 * step)) / (2 * delta)
    } else {
      (fn(par + step) - fn(par - step)) / (2 * delta)
    }
  }

  hessian <- matrix(NA_real_, n, n)
  for (j in free_idx) {
    step <- rep(0, n)
    step[j] <- delta
    hessian[j, j] <- (fn(par + step) - 2 * f0 + fn(par - step)) / delta^2
  }
  if (length(free_idx) > 1L) {
    for (a in 2:length(free_idx)) {
      for (b in seq_len(a - 1L)) {
        i <- free_idx[a]
        j <- free_idx[b]
        si <- sj <- rep(0, n)
        si[i] <- delta
        sj[j] <- delta
        hessian[i, j] <- hessian[j, i] <-
          (fn(par + si + sj) - fn(par + si - sj) -
            fn(par - si + sj) + fn(par - si - sj)) / (4 * delta^2)
      }
    }
  }

  h_free <- hessian[free_idx, free_idx, drop = FALSE]
  g_free <- gradient[free_idx]
  eigenvalues <- tryCatch(
    eigen(h_free, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) rep(NA_real_, length(free_idx))
  )
  scaled <- tryCatch(solve(chol(h_free), g_free),
    error = function(e) rep(NA_real_, length(free_idx))
  )
  max_raw <- if (length(g_free)) max(abs(g_free)) else 0
  max_scaled <- if (length(scaled) && all(is.finite(scaled))) {
    max(abs(scaled))
  } else {
    NA_real_
  }
  max_min <- if (length(g_free) && length(scaled) && all(is.finite(scaled))) {
    max(pmin(abs(g_free), abs(scaled)))
  } else {
    NA_real_
  }

  kkt_margin <- c(gradient[active_lower], -gradient[active_upper])
  kkt_min <- if (length(kkt_margin)) min(kkt_margin) else Inf
  kkt_ok <- all(is.finite(kkt_margin)) && kkt_min >= -gradient_tol
  h_min <- if (length(eigenvalues)) min(eigenvalues) else Inf
  h_max <- if (length(eigenvalues)) max(eigenvalues) else Inf
  h_ratio <- if (is.finite(h_min) && is.finite(h_max) && h_max > 0) {
    h_min / h_max
  } else {
    NA_real_
  }
  hessian_ok <- all(is.finite(eigenvalues)) &&
    h_min > hessian_eigen_tol && is.finite(h_ratio) &&
    h_ratio >= hessian_ratio_tol
  gradient_ok <- is.finite(max_min) && max_min <= gradient_tol

  data.frame(
    objective = f0, n_parameters = n, n_free = length(free_idx),
    n_active_lower = sum(active_lower), n_active_upper = sum(active_upper),
    max_abs_gradient = max_raw, max_abs_scaled_gradient = max_scaled,
    max_min_gradient = max_min, gradient_ok = gradient_ok,
    kkt_min_margin = kkt_min, kkt_ok = kkt_ok,
    hessian_min_eigen = h_min, hessian_max_eigen = h_max,
    hessian_eigen_ratio = h_ratio, hessian_ok = hessian_ok,
    numerical_ok = gradient_ok && kkt_ok && hessian_ok,
    stringsAsFactors = FALSE
  )
}

ilhte_mermod_diagnostics <- function(model, ...) {
  stopifnot(inherits(model, "merMod"))
  theta <- lme4::getME(model, "theta")
  lower_theta <- lme4::getME(model, "lower")
  par <- c(theta, lme4::fixef(model))
  lower <- c(lower_theta, rep(-Inf, length(lme4::fixef(model))))
  objective <- lme4::getME(model, "devfun")
  nd <- ilhte_objective_diagnostics(objective, par, lower = lower, ...)

  opt_code <- model@optinfo$conv$opt
  optimizer_ok <- length(opt_code) > 0L && all(is.finite(opt_code)) &&
    all(opt_code == 0)
  messages <- as.character(model@optinfo$conv$lme4$messages)
  boundary_message <- grepl("boundary \\(singular\\)", messages,
    ignore.case = TRUE
  )
  nonboundary <- messages[!boundary_message]
  vc <- as.data.frame(lme4::VarCorr(model))
  variance_rows <- is.na(vc$var2)
  slope_row <- vc$grp == "item" & vc$var1 == "treat" & variance_rows
  corr_row <- vc$grp == "item" & !is.na(vc$var2) &
    ((vc$var1 == "(Intercept)" & vc$var2 == "treat") |
      (vc$var2 == "(Intercept)" & vc$var1 == "treat"))
  slope_sd <- if (any(slope_row)) vc$sdcor[which(slope_row)[1]] else NA_real_
  slope_variance <- slope_sd^2
  item_correlation <- if (any(corr_row)) vc$sdcor[which(corr_row)[1]] else NA_real_
  finite_bounds <- is.finite(lower_theta)
  boundary_cholesky <- any(finite_bounds &
    abs(theta - lower_theta) <= 1e-7)
  boundary_variance <- any(variance_rows & is.finite(vc$vcov) & vc$vcov <= 1e-8)
  boundary_correlation <- is.finite(item_correlation) &&
    abs(item_correlation) >= 1 - 1e-6

  cbind(data.frame(
    optimizer_code = if (length(opt_code)) paste(opt_code, collapse = ";") else NA,
    optimizer_ok = optimizer_ok,
    lme4_messages = paste(messages, collapse = " | "),
    nonboundary_message_count = length(nonboundary),
    nonboundary_messages = paste(nonboundary, collapse = " | "),
    boundary_cholesky = boundary_cholesky,
    boundary_variance = boundary_variance,
    boundary_slope_variance = is.finite(slope_variance) && slope_variance <= 1e-8,
    boundary_correlation = boundary_correlation,
    singular = lme4::isSingular(model, tol = 1e-4),
    random_slope_sd = slope_sd,
    random_slope_variance = slope_variance,
    item_intercept_slope_correlation = item_correlation,
    stringsAsFactors = FALSE
  ), nd)
}

ilhte_replicate <- function(n_sub, n_items, te_mean, te_sd, disc = 1.7,
                            seed = 1L) {
  d <- ilhte_data(n_sub, n_items, te_mean, te_sd, disc, seed)
  r <- ilhte_fit(d)
  cbind(data.frame(
    n_sub = n_sub, n_items = n_items, te_mean = te_mean,
    te_sd = te_sd, disc = disc, seed = seed,
    stringsAsFactors = FALSE
  ), r)
}

# Find the common discrimination needed for this item bank and target.
ilhte_disc_for_rho <- function(b, target, theta) {
  f <- function(disc) {
    rho_three(1, theta, b, rep(disc, length(b)),
      theta_var = 1
    )$rho_tilde - target
  }
  if (f(0.05) > 0 || f(20) < 0) {
    return(NA_real_)
  }
  stats::uniroot(f, c(0.05, 20), tol = 1e-8)$root
}

# Calibrate before drawing responses. Under theta_fixed, scale the logit
# treatment effects with discrimination to hold their latent-coordinate size fixed.
ilhte_replicate_targeted <- function(n_sub, n_items, target_rho, te_mean, te_sd,
                                     te_scale = c("theta_fixed", "logit_fixed"),
                                     theta_quad = NULL, seed = 1L,
                                     disc_reference = 1.7) {
  te_scale <- match.arg(te_scale)
  set.seed(seed)
  b <- stats::rnorm(n_items)
  if (is.null(theta_quad)) theta_quad <- stats::rnorm(20000L)

  disc <- ilhte_disc_for_rho(b, target_rho, theta_quad)
  if (is.na(disc)) {
    return(cbind(
      data.frame(
        n_sub = n_sub, n_items = n_items,
        target_rho = target_rho, te_mean = te_mean,
        te_sd = te_sd, te_scale = te_scale,
        disc = NA_real_, realized_rho = NA_real_,
        seed = seed, stringsAsFactors = FALSE
      ),
      data.frame(
        method = c("Constant", "Varying"),
        estimate = NA_real_, std.error = NA_real_,
        statistic = NA_real_, p.value = NA_real_,
        converged = FALSE, sig_het = NA_real_,
        stringsAsFactors = FALSE
      )
    ))
  }
  realized <- rho_three(1, theta_quad, b, rep(disc, n_items),
    theta_var = 1
  )$rho_tilde

  k <- if (te_scale == "theta_fixed") disc / disc_reference else 1
  d <- ilhte_data_with_bank(n_sub, b, te_mean * k, te_sd * k, disc, seed)
  r <- ilhte_fit(d)
  cbind(data.frame(
    n_sub = n_sub, n_items = n_items, target_rho = target_rho,
    te_mean = te_mean, te_sd = te_sd, te_scale = te_scale,
    disc = disc, realized_rho = realized, seed = seed,
    stringsAsFactors = FALSE
  ), r)
}

ilhte_data_with_bank <- function(n_sub, b, te_mean, te_sd, disc, seed) {
  n_items <- length(b)

  set.seed(seed + 500000L)
  theta <- stats::rnorm(n_sub)
  treat <- rep(0L, n_sub)
  treat[sample.int(n_sub, floor(n_sub / 2))] <- 1L
  set.seed(seed + 700000L)
  tau <- stats::rnorm(n_items, te_mean, te_sd)
  u <- matrix(stats::runif(n_sub * n_items), nrow = n_sub)

  lin <- outer(theta, b, "-") * disc + outer(treat, tau, "*")
  y <- 1L * (u < stats::plogis(lin))
  data.frame(
    person = rep(seq_len(n_sub), times = n_items),
    item = rep(seq_len(n_items), each = n_sub),
    treat = rep(treat, times = n_items),
    correct = as.vector(y), stringsAsFactors = FALSE
  )
}

.ilhte_key_piece <- function(x) {
  x <- enc2utf8(as.character(x))
  paste0(nchar(x, type = "bytes"), ":", x)
}

.ilhte_identity_key <- function(x, columns) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) {
    stop("retry output is missing identity columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyNA(x[columns])) {
    stop("retry identity columns must not contain missing values",
      call. = FALSE
    )
  }
  for (column in columns) {
    if (any(!nzchar(trimws(as.character(x[[column]]))))) {
      stop("retry identity column contains a blank value: ", column,
        call. = FALSE
      )
    }
  }
  vapply(
    seq_len(nrow(x)), function(i) {
      paste(vapply(columns, function(column) {
        .ilhte_key_piece(x[[column]][i])
      }, character(1)), collapse = "|")
    },
    character(1)
  )
}

.ilhte_row_signature <- function(x) {
  vapply(seq_len(nrow(x)), function(i) {
    paste(vapply(names(x), function(column) {
      z <- x[[column]][i]
      value <- if (is.na(z)) "<NA>" else as.character(z)
      paste0(.ilhte_key_piece(column), "=", .ilhte_key_piece(value))
    }, character(1)), collapse = "|")
  }, character(1))
}

assert_retry_output_identity <- function(
  x, bank_columns, attempt_column = "attempt",
  expected_attempts = NULL
) {
  if (!is.data.frame(x) || !nrow(x)) {
    stop("retry output must be a non-empty data frame", call. = FALSE)
  }
  identity_columns <- c(bank_columns, attempt_column)
  identity <- .ilhte_identity_key(x, identity_columns)
  if (anyDuplicated(identity)) {
    stop("retry output has more than one row for a bank/attempt identity",
      call. = FALSE
    )
  }
  bank <- .ilhte_identity_key(x, bank_columns)
  if (!is.null(expected_attempts)) {
    if (length(expected_attempts) == 1L && is.numeric(expected_attempts)) {
      expected <- as.integer(expected_attempts)
      if (is.na(expected) || expected < 1L ||
        any(tabulate(match(bank, unique(bank))) != expected)) {
        stop("retry output does not contain the expected attempt count per bank",
          call. = FALSE
        )
      }
    } else {
      expected <- sort(as.character(expected_attempts), method = "radix")
      actual <- split(as.character(x[[attempt_column]]), bank)
      if (any(vapply(actual, function(z) {
        !identical(sort(z, method = "radix"), expected)
      }, logical(1)))) {
        stop("retry output does not contain the declared attempt set per bank",
          call. = FALSE
        )
      }
    }
  }
  invisible(TRUE)
}

deduplicate_retry_attempts <- function(
  x, bank_columns, attempt_column = "attempt",
  expected_attempts = NULL
) {
  if (!is.data.frame(x) || !nrow(x)) {
    stop("retry output must be a non-empty data frame", call. = FALSE)
  }
  identity <- .ilhte_identity_key(x, c(bank_columns, attempt_column))
  signature <- .ilhte_row_signature(x)
  groups <- split(seq_len(nrow(x)), identity)
  conflicting <- vapply(groups, function(i) {
    length(unique(signature[i])) != 1L
  }, logical(1))
  if (any(conflicting)) {
    stop("conflicting rows share a bank/attempt retry identity",
      call. = FALSE
    )
  }
  keep <- !duplicated(identity)
  out <- x[keep, , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "duplicate_rows_removed") <- sum(!keep)
  assert_retry_output_identity(
    out, bank_columns, attempt_column,
    expected_attempts
  )
  out
}

deduplicate_bank_retry_counts <- function(
  x, bank_columns, retry_column = "bank_retries"
) {
  if (!is.data.frame(x) || !nrow(x)) {
    stop("bank retry output must be a non-empty data frame", call. = FALSE)
  }
  if (!retry_column %in% names(x)) {
    stop("bank retry output is missing ", retry_column, call. = FALSE)
  }
  bank <- .ilhte_identity_key(x, bank_columns)
  retry <- x[[retry_column]]
  if (!is.numeric(retry) || anyNA(retry) || any(!is.finite(retry)) ||
    any(retry < 0 | retry != floor(retry))) {
    stop("bank retry counts must be finite non-negative integers",
      call. = FALSE
    )
  }
  groups <- split(seq_len(nrow(x)), bank)
  if (any(vapply(
    groups, function(i) length(unique(retry[i])) != 1L,
    logical(1)
  ))) {
    stop("conflicting retry counts share a bank identity", call. = FALSE)
  }

  first <- which(!duplicated(bank))
  out <- x[first, c(bank_columns, retry_column), drop = FALSE]
  rownames(out) <- NULL
  counts <- table(factor(bank, levels = bank[first]))
  out$rows_collapsed <- as.integer(counts)
  if (anyDuplicated(.ilhte_identity_key(out, bank_columns))) {
    stop("bank retry deduplication did not produce one row per bank",
      call. = FALSE
    )
  }
  out
}
