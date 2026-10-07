# Quadrature rules and deterministic checks of information-index behavior.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

theory_grid <- function(shape = "normal", n_node = 12001L, range = 8,
                        df = 5) {
  theta <- seq(-range, range, length.out = n_node)
  dens <- switch(shape,
    normal = stats::dnorm(theta),
    skew_pos = {
      k <- 4
      s <- sqrt(k)
      ifelse(theta * s + k > 0,
        stats::dgamma(theta * s + k, shape = k, rate = 1) * s,
        0
      )
    },
    bimodal = 0.5 * stats::dnorm(theta, -1, 0.6) +
      0.5 * stats::dnorm(theta, 1, 0.6),
    heavy_tail = stats::dt(theta, df = df),
    stop("unknown shape: ", shape)
  )

  w <- dens / sum(dens)
  mu <- sum(w * theta)
  v <- sum(w * (theta - mu)^2)

  theta <- (theta - mu) / sqrt(v)

  excluded <- switch(shape,
    normal     = 2 * stats::pnorm(-range),
    skew_pos   = stats::pgamma(-range * sqrt(4) + 4, shape = 4, rate = 1),
    bimodal    = stats::pnorm(-range, -1, 0.6) + stats::pnorm(-range, 1, 0.6),
    heavy_tail = 2 * stats::pt(-range, df = df)
  )
  kurt <- sum(w * (theta - sum(w * theta))^4)

  list(
    theta = theta, weights = w, var = 1, shape = shape, n_node = n_node,
    range = range, tail_mass_excluded = excluded, kurtosis = kurt
  )
}

rho_curve3 <- function(cs, grid, beta, lambda0, guessing = NULL,
                       tail_integrable = TRUE) {
  out <- lapply(cs, function(cc) {
    r <- rho_three(cc, grid$theta, beta, lambda0,
      guessing = guessing,
      theta_var = 1, weights = grid$weights,
      tail_integrable = tail_integrable,
      wbar_policy = "value"
    )
    cbind(data.frame(c = cc), r)
  })
  do.call(rbind, out)
}

count_roots <- function(cs, values, level) {
  d <- values - level
  d <- d[is.finite(d)]
  if (length(d) < 2L) {
    return(NA_integer_)
  }
  sum(d[-1] * d[-length(d)] < 0)
}

monotone_report <- function(values, tol = 1e-10) {
  d <- diff(values)
  d <- d[is.finite(d)]
  list(
    monotone = all(d > -tol),
    n_decrease = sum(d < -tol),
    max_decrease = if (any(d < 0)) -min(d) else 0
  )
}

phi_root <- function() {
  stats::uniroot(function(x) 2 - x * tanh(x / 2), c(1, 6), tol = .Machine$double.eps^0.75)$root
}

theory_grid_adaptive <- function(shape = "normal", beta, c, lambda0 = 1,
                                 n_global = 20001L, n_local = 4001L,
                                 range = 8, df = 5) {
  raw_density <- function(x) {
    switch(shape,
      normal = stats::dnorm(x),
      skew_pos = {
        k <- 4
        s <- sqrt(k)
        ifelse(x * s + k > 0,
          stats::dgamma(x * s + k, shape = k, rate = 1) * s, 0
        )
      },
      bimodal = 0.5 * stats::dnorm(x, -1, 0.6) + 0.5 * stats::dnorm(x, 1, 0.6),
      heavy_tail = stats::dt(x, df = df),
      stop("unknown shape: ", shape)
    )
  }

  t0 <- seq(-range, range, length.out = n_global)
  w0 <- raw_density(t0)
  w0 <- w0 / sum(w0)
  mu_raw <- sum(w0 * t0)
  s_raw <- sqrt(sum(w0 * (t0 - mu_raw)^2))
  beta_raw <- mu_raw + s_raw * beta

  lam <- max(lambda0)
  half_std <- max(6 * log1p(c * lam) / (c * lam), 1e-6)
  half_raw <- half_std * s_raw

  theta <- sort(unique(c(
    t0,
    unlist(lapply(beta_raw, function(b) {
      seq(b - half_raw, b + half_raw, length.out = n_local)
    }))
  )))
  theta <- theta[theta >= -range & theta <= range]

  dens <- raw_density(theta)
  n <- length(theta)
  dl <- c(theta[2] - theta[1], diff(theta))
  dr <- c(diff(theta), theta[n] - theta[n - 1])
  w <- dens * (dl + dr) / 2
  w <- w / sum(w)

  mu <- sum(w * theta)
  v <- sum(w * (theta - mu)^2)
  list(
    theta = (theta - mu) / sqrt(v), weights = w, var = 1, shape = shape,
    n_node = n, half_width = half_std, range = range, raw_scale = s_raw
  )
}

theory_density_std <- function(shape, df = 5) {
  s <- switch(shape,
    normal = 1,
    skew_pos = 1,
    bimodal = sqrt(1 + 0.36),
    heavy_tail = sqrt(df / (df - 2)),
    stop("unknown shape: ", shape)
  )
  function(x) {
    y <- s * x
    s * switch(shape,
      normal = stats::dnorm(y),
      skew_pos = {
        k <- 4
        sq <- sqrt(k)
        ifelse(y * sq + k > 0,
          stats::dgamma(y * sq + k, shape = k, rate = 1) * sq, 0
        )
      },
      bimodal = 0.5 * stats::dnorm(y, -1, 0.6) + 0.5 * stats::dnorm(y, 1, 0.6),
      heavy_tail = stats::dt(y, df = df)
    )
  }
}

logistic_variance_kernel <- function(x) {
  a <- exp(-abs(x))
  a / (1 + a)^2
}

normal_mixture_2pl_mean_information <- function(
  c, weight, mean, sd, beta = 0, lambda0 = 1, rel.tol = 1e-10
) {
  stopifnot(
    is.numeric(c), length(c) >= 1L, all(is.finite(c)), all(c > 0),
    is.numeric(weight), is.numeric(mean), is.numeric(sd),
    length(weight) == length(mean), length(mean) == length(sd),
    all(is.finite(weight)), all(weight >= 0), sum(weight) > 0,
    all(is.finite(mean)), all(is.finite(sd)), all(sd > 0),
    length(beta) == 1L, is.finite(beta),
    length(lambda0) == 1L, is.finite(lambda0), lambda0 > 0
  )
  weight <- weight / sum(weight)

  vapply(c, function(cc) {
    a <- cc * lambda0
    sum(weight * vapply(seq_along(weight), function(k) {
      stats::integrate(
        function(z) {
          x <- a * (mean[k] + sd[k] * z - beta)
          a^2 * logistic_variance_kernel(x) * stats::dnorm(z)
        },
        -Inf, Inf,
        subdivisions = 3000L, rel.tol = rel.tol
      )$value
    }, numeric(1)))
  }, numeric(1))
}

smooth_three_normal_counterexample <- function(level = 0.02) {
  weight <- c(0.9998, 0.0001, 0.0001)
  mean <- c(3, 0, -100)
  sd <- rep(0.15, 3)
  c_grid <- c(0.1, 0.3, 0.5, 0.8, 1, 1.5, 2, 4, 10, 50, 100)
  brackets <- rbind(c(0.3, 0.5), c(1.5, 2), c(50, 100))
  jbar <- function(z) {
    normal_mixture_2pl_mean_information(
      z,
      weight = weight, mean = mean, sd = sd
    )
  }
  roots <- apply(brackets, 1L, function(b) {
    stats::uniroot(function(z) jbar(z) - level, b, tol = 1e-10)$root
  })
  latent_mean <- sum(weight * mean)
  latent_variance <- sum(weight * (sd^2 + mean^2)) - latent_mean^2
  density_at_beta <- sum(weight * stats::dnorm(0, mean, sd))

  list(
    curve = data.frame(c = c_grid, Jbar = jbar(c_grid)),
    roots = as.numeric(roots),
    brackets = brackets,
    level = level,
    latent_mean = latent_mean,
    latent_variance = latent_variance,
    reliability_target = latent_variance * level /
      (1 + latent_variance * level),
    density_at_beta = density_at_beta,
    density_upper_bound = sum(weight / (sd * sqrt(2 * pi)))
  )
}

one_item_normal_wbar_exact <- function(c) {
  stopifnot(is.numeric(c), length(c) >= 1L, all(is.finite(c)), all(c > 0))
  log1pexp <- function(x) {
    ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
  }
  z <- c^2 / 2
  log_msem <- log(2) + log1pexp(z) - 2 * log(c)
  log_wbar <- -log1pexp(log_msem)
  log_asymptotic <- 2 * log(c) - log(2) - z
  data.frame(
    c = c,
    log_msem = log_msem,
    log_wbar = log_wbar,
    log_wbar_over_c = log_wbar / c,
    log_asymptotic = log_asymptotic,
    asymptotic_ratio = exp(log_wbar - log_asymptotic)
  )
}

uniform_desert_reliability <- function(c) {
  stopifnot(is.numeric(c), length(c) >= 1L, all(is.finite(c)), all(c > 0))

  a <- exp(-c)
  jbar <- c * a * (1 - a) / ((1 + a^2) * (1 + a))
  scaled <- jbar / 12
  data.frame(
    c = c,
    Jbar = jbar,
    asymptotic_ratio = jbar / (c * exp(-c)),
    rho_tilde = scaled / (1 + scaled)
  )
}

reciprocal_information_tail_lower_log <- function(
  R, c = 1, lambda = 1, beta = 0, c0 = 2, nu = 2
) {
  stopifnot(
    is.numeric(R), length(R) >= 1L, all(is.finite(R)), all(R > 0),
    length(c) == 1L, is.finite(c), c > 0,
    length(lambda) == 1L, is.finite(lambda), lambda > 0,
    length(beta) == 1L, is.finite(beta),
    length(c0) == 1L, is.finite(c0), c0 > 0,
    length(nu) == 1L, is.finite(nu), nu > 0
  )
  a <- c * lambda
  log(c0) - 2 * log(a) + a * (R - beta) - (nu + 1) * log(R + 1)
}

guessing_spike_coefficient <- function(g) {
  stopifnot(
    is.numeric(g), length(g) >= 1L, all(is.finite(g)),
    all(g >= 0), all(g < 1)
  )
  x <- 1 - g
  out <- numeric(length(g))
  out[g == 0] <- 1
  regular <- g > 0 & x >= 1e-5
  out[regular] <- 1 + g[regular] / x[regular] * log(g[regular])
  near_one <- g > 0 & x < 1e-5

  xx <- x[near_one]
  out[near_one] <- xx / 2 + xx^2 / 6 + xx^3 / 12 + xx^4 / 20
  out
}
