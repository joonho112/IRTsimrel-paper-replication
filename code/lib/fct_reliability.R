# Stable test-information calculations and the three information indices.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

# Compute log test information without subtracting probabilities near one.
# Work in chunks to bound memory when a study uses many integration nodes.
log_test_information <- function(c, theta, beta, lambda0, guessing = NULL,
                                 chunk_size = NULL) {
  stopifnot(is.numeric(c), length(c) == 1L, is.finite(c), c > 0)
  stopifnot(is.numeric(theta), length(theta) >= 1L, all(is.finite(theta)))
  stopifnot(is.numeric(beta), length(beta) >= 1L, all(is.finite(beta)))

  n_item <- length(beta)
  if (length(lambda0) == 1L) lambda0 <- rep(lambda0, n_item)
  stopifnot(length(lambda0) == n_item, all(is.finite(lambda0)), all(lambda0 > 0))

  if (is.null(guessing)) guessing <- 0
  if (length(guessing) == 1L) guessing <- rep(guessing, n_item)
  stopifnot(
    length(guessing) == n_item, all(is.finite(guessing)),
    all(guessing >= 0), all(guessing < 1)
  )

  lambda <- c * lambda0
  log_lambda2 <- 2 * log(lambda)
  log_1mg <- log1p(-guessing)
  any_g <- any(guessing > 0)

  if (is.null(chunk_size)) {
    chunk_size <- max(1L, min(length(theta), as.integer(4e6 / n_item)))
  }

  out <- numeric(length(theta))
  starts <- seq.int(1L, length(theta), by = chunk_size)
  for (s in starts) {
    idx <- s:min(s + chunk_size - 1L, length(theta))
    eta <- outer(theta[idx], beta, "-") *
      rep(lambda, each = length(idx))
    log_p_star <- stats::plogis(eta, log.p = TRUE)
    log_q_star <- stats::plogis(-eta, log.p = TRUE)

    if (!any_g) {
      log_j <- rep(log_lambda2, each = length(idx)) + log_p_star + log_q_star
    } else {
      p <- rep(guessing, each = length(idx)) +
        rep(1 - guessing, each = length(idx)) * exp(log_p_star)
      log_p <- ifelse(rep(guessing, each = length(idx)) > 0,
        log(p), log_p_star
      )
      log_j <- rep(log_lambda2, each = length(idx)) +
        rep(log_1mg, each = length(idx)) +
        2 * log_p_star + log_q_star - log_p
    }
    out[idx] <- .log_row_sum_exp(log_j)
  }
  out
}

.log_row_sum_exp <- function(m) {
  pivot <- apply(m, 1L, max)
  pivot[!is.finite(pivot)] <- 0
  res <- pivot + log(rowSums(exp(m - pivot)))
  res[is.nan(res)] <- -Inf
  res
}

.log_mean_exp <- function(x, log_w) {
  keep <- is.finite(log_w)
  x <- x[keep]
  log_w <- log_w[keep]
  if (any(is.infinite(x) & x > 0)) {
    return(Inf)
  }
  z <- x + log_w
  fin <- z[is.finite(z)]
  if (!length(fin)) {
    return(-Inf)
  }
  pivot <- max(fin)
  pivot + log(sum(exp(z - pivot)))
}

# Use the supplied variance, or the variance represented by the nodes.
# Equal Monte Carlo weights use sample variance; explicit weights use a population moment.
resolve_theta_var <- function(theta, theta_var = NULL, weights = NULL) {
  if (!is.null(theta_var)) {
    stopifnot(
      is.numeric(theta_var), length(theta_var) == 1L,
      is.finite(theta_var), theta_var > 0
    )
    return(structure(as.numeric(theta_var), source = "supplied"))
  }
  if (is.null(weights)) {
    v <- stats::var(theta)
    stopifnot(is.finite(v), v > 0)
    return(structure(v, source = "sample_variance"))
  }
  w <- weights / sum(weights)
  mu <- sum(w * theta)
  v <- sum(w * (theta - mu)^2)
  stopifnot(is.finite(v), v > 0)
  structure(v, source = "weighted_population_variance")
}

rho_at_expected_information <- function(mean_information, theta_var) {
  stopifnot(
    is.numeric(mean_information), length(mean_information) == 1L,
    is.finite(mean_information), mean_information > 0,
    is.numeric(theta_var), length(theta_var) == 1L,
    is.finite(theta_var), theta_var > 0
  )
  stats::plogis(log(theta_var) + log(mean_information))
}

# Keep the primary transform of mean information separate from the
# descriptive mean of the form-specific indices. Jensen's inequality separates them.
item_superpopulation_reliability <- function(form_mean_information, theta_var,
                                             form_reliability = NULL) {
  stopifnot(
    is.numeric(form_mean_information),
    length(form_mean_information) >= 1L,
    all(is.finite(form_mean_information)),
    all(form_mean_information > 0),
    is.numeric(theta_var), length(theta_var) == 1L,
    is.finite(theta_var), theta_var > 0
  )
  expected_form_reliability <- stats::plogis(
    log(theta_var) + log(form_mean_information)
  )
  if (is.null(form_reliability)) {
    form_reliability <- expected_form_reliability
  } else {
    stopifnot(
      is.numeric(form_reliability),
      length(form_reliability) == length(form_mean_information),
      all(is.finite(form_reliability)),
      all(form_reliability >= 0), all(form_reliability <= 1)
    )
    if (max(abs(form_reliability - expected_form_reliability)) > 1e-10) {
      stop("Form reliabilities do not correspond to the supplied form mean ",
        "information and theta variance.",
        call. = FALSE
      )
    }
  }

  expected_information <- mean(form_mean_information)
  primary <- rho_at_expected_information(expected_information, theta_var)
  descriptive <- mean(form_reliability)
  data.frame(
    expected_information = expected_information,
    rho_at_expected_information = primary,
    mean_random_form_reliability = descriptive,
    jensen_gap = primary - descriptive,
    n_forms = length(form_mean_information),
    stringsAsFactors = FALSE
  )
}

# Evaluate all three indices on the same form and integration rule.
# A finite tail truncation is stored separately when the population reciprocal
# information expectation does not exist.
rho_three <- function(c, theta, beta, lambda0, guessing = NULL,
                      theta_var = NULL, weights = NULL,
                      tail_integrable = TRUE,
                      wbar_policy = c("auto", "value")) {
  wbar_policy <- match.arg(wbar_policy)
  n <- length(theta)
  log_j <- log_test_information(c, theta, beta, lambda0, guessing)

  if (is.null(weights)) {
    log_w <- rep(-log(n), n)
  } else {
    stopifnot(
      is.numeric(weights), length(weights) == n,
      all(is.finite(weights)), all(weights >= 0), sum(weights) > 0
    )
    pos <- weights > 0
    lw <- rep(-Inf, n)
    lw[pos] <- log(weights[pos])
    pivot <- max(lw[pos])
    log_w <- lw - (pivot + log(sum(exp(lw[pos] - pivot))))
  }

  sv <- resolve_theta_var(theta, theta_var, weights)
  log_sv <- log(as.numeric(sv))

  log_jbar <- .log_mean_exp(log_j, log_w)
  log_msem <- .log_mean_exp(-log_j, log_w)

  rho_tilde <- stats::plogis(log_sv + log_jbar)

  active <- is.finite(log_w)
  rho_psd <- sum(exp(log_w[active]) * stats::plogis(log_sv + log_j[active]))
  w_bar_raw <- stats::plogis(log_sv - log_msem)

  zero_info <- any(is.infinite(log_j[active]) & log_j[active] < 0)

  log_terms <- (-log_j[active]) + log_w[active]
  msem_top_share <- if (any(is.finite(log_terms))) {
    pivot <- max(log_terms[is.finite(log_terms)])
    exp(pivot - .log_mean_exp(-log_j, log_w))
  } else {
    NA_real_
  }

  status <- if (zero_info) {
    "zero_information"
  } else if (!tail_integrable) {
    "wbar_undefined"
  } else {
    "ok"
  }
  w_bar <- if (wbar_policy == "auto" && status != "ok") NA_real_ else w_bar_raw

  data.frame(
    rho_tilde = rho_tilde,
    rho_psd = rho_psd,
    w_bar = w_bar,
    w_bar_truncated = w_bar_raw,
    mean_info = exp(log_jbar),
    msem = exp(log_msem),
    msem_top_share = msem_top_share,
    theta_var = as.numeric(sv),
    theta_var_source = attr(sv, "source"),
    n_nodes = n,
    status = status,
    stringsAsFactors = FALSE
  )
}

rho_psd <- function(c, theta, beta, lambda0, guessing = NULL,
                    theta_var = NULL, weights = NULL) {
  rho_three(c, theta, beta, lambda0, guessing, theta_var, weights)$rho_psd
}
