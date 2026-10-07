# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

# Section 4: choose a target, calibrate a form, then evaluate fresh draws.
source("code/00_setup.R")
use_irtsimrel("0.3.1")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
fit <- IRTsimrel::eqc_calibrate(
  target_rho = .85, n_items = 20L, model = "2pl",
  latent_shape = "normal", item_source = "parametric",
  item_params = list(difficulty_params = list(mu = 0, sigma = 1)),
  reliability_metric = "info", M = 20000L, c_bounds = c(.1, 10),
  tol = 1e-6, root_policy = "lowest_increasing", seed = 20260826L
)
# The independent evaluation represents variance one. The solver uses the
# variance represented by its own integration nodes, recorded below.
set.seed(20260827L)
theta <- rnorm(200000L)
theta <- (theta - mean(theta)) / sd(theta)
fresh <- IRTsimrel::compute_rho_tilde(1, theta, fit$beta_vec,
  fit$lambda_scaled,
  theta_var = 1
)
answer <- data.frame(
  target = .85, multiplier = fit$c_star,
  solver_index = fit$achieved_rho, solver_variance = var(as.numeric(fit$theta_quad)),
  fresh_index = fresh, status = fit$misc$calibration_status
)
write_result(answer, repo_path("output", "example", "calibration.csv"))
write_result(data.frame(
  item = seq_along(fit$beta_vec), difficulty = fit$beta_vec,
  discrimination = fit$lambda_scaled
), repo_path("output", "example", "item_form.csv"))
reference <- read_result("software", "normal_calibration.csv")
stopifnot(
  abs(fit$c_star - reference$multiplier[1]) < 1e-8,
  abs(fit$achieved_rho - .85) < 1e-6
)
record_session(repo_path("output", "example", "session-info.txt"))
print(answer)
