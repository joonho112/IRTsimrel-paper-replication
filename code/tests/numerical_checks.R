# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

# These checks use analytic identities and deliberately invalid inputs.
# They are sourced by code/03_verify.R after the reliability functions.
beta <- c(-1, 0, .6, 1.2)
lambda <- c(.7, 1, 1.3, 1.8)
theta <- seq(-4, 4, length.out = 1001)
w <- dnorm(theta)
w <- w / sum(w)
r <- rho_three(1.2, theta, beta, lambda, theta_var = 1, weights = w)
check(
  "ordering of the three finite-grid indices",
  as.integer(r$w_bar <= r$rho_psd && r$rho_psd <= r$rho_tilde), 1
)
# Changing latent coordinates should preserve each index when variance and
# discrimination are transformed together.
s <- rho_three(1.2, 3 * theta + 2, 3 * beta + 2, lambda / 3,
  theta_var = 9, weights = w
)
for (name in c("rho_tilde", "rho_psd", "w_bar")) {
  check(paste("affine invariance", name), s[[name]], r[[name]], 1e-12)
}
# For a single 2PL item, information is a^2 p(1-p).
a <- 1.4
th <- seq(-3, 3, length.out = 41)
observed <- exp(log_test_information(1, th, .2, a))
p <- plogis(a * (th - .2))
check("single-item information formula", max(abs(observed - a^2 * p * (1 - p))), 0, 1e-14)
extreme <- log_test_information(1, c(-1000, 1000), c(-1, 1), c(1, 2))
check("log information remains finite in extreme tails", sum(is.finite(extreme)), 2)
tail <- rho_three(1, theta, beta, lambda,
  theta_var = 1,
  weights = w, tail_integrable = FALSE
)
check(
  "divergent-tail integral is not reported as a finite population index",
  as.integer(is.na(tail$w_bar) && is.finite(tail$w_bar_truncated)), 1
)
z <- item_superpopulation_reliability(c(.1, 1, 10), 1)
check(
  "transform of mean information differs from mean transformed information",
  as.integer(z$jensen_gap > 0), 1
)
fails <- function(expr) {
  inherits(tryCatch(
    {
      force(expr)
      NULL
    },
    error = identity
  ), "error")
}
check("invalid discrimination is rejected", as.integer(fails(
  log_test_information(1, theta, beta, c(1, 0, 1, 1))
)), 1)
check("inconsistent form indices are rejected", as.integer(fails(
  item_superpopulation_reliability(c(.1, 1, 10), 1, c(.9, .9, .9))
)), 1)
