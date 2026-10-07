# Item-generation recipes and their induced information distributions.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

# Translate the documented recipe into item-generation arguments.
recipe_to_args <- function(recipe, pool = c("parametric", "irw"),
                           pool_data = NULL) {
  pool <- match.arg(pool)
  stopifnot(nrow(recipe) == 1L)

  args <- list(
    n_items = as.integer(recipe$n_items),
    model = recipe$model,
    center_difficulties = TRUE
  )

  if (recipe$model == "rasch" && !is.na(recipe$lambda_fixed) &&
    recipe$lambda_fixed != 1) {
    args$scale <- as.numeric(recipe$lambda_fixed)
  }

  if (pool == "irw") {
    if (is.null(pool_data)) stop("pool = 'irw' needs pool_data.", call. = FALSE)
    args$source <- "irw"
    args$difficulty_params <- list(pool = pool_data)
  } else {
    args$source <- "parametric"
    args$difficulty_params <- list(
      distribution = recipe$beta_dist,
      mu = as.numeric(recipe$beta_mu),
      sigma = as.numeric(recipe$beta_sigma)
    )
  }

  if (recipe$model != "rasch") {
    args$discrimination_params <- switch(recipe$lambda_spec,
      lognormal = list(
        mu_log = as.numeric(recipe$lambda_mu_log),
        sigma_log = as.numeric(recipe$lambda_sigma_log),
        rho = if (is.na(recipe$ab_rho)) 0 else as.numeric(recipe$ab_rho)
      ),
      fixed = list(
        mu_log = log(as.numeric(recipe$lambda_fixed)),
        sigma_log = 1e-8, rho = 0
      ),
      uniform = list(mu_log = 0.16628, sigma_log = 0.33707, rho = 0),
      stop("unknown lambda_spec: ", recipe$lambda_spec)
    )
    args$method <- if (!is.na(recipe$ab_rho) && recipe$ab_rho != 0) {
      "copula"
    } else {
      "independent"
    }
  }

  if (recipe$model == "3pl") {
    args$guessing_params <- switch(recipe$guessing_dist,
      fixed   = list(distribution = "fixed", value = as.numeric(recipe$guessing_value)),
      beta    = list(distribution = "beta", shape1 = 5, shape2 = 17),
      uniform = list(distribution = "uniform", min = 0.10, max = 0.30),
      stop("unknown guessing_dist: ", recipe$guessing_dist)
    )
  }
  args
}

# Draw independent item forms from one recipe, keeping the integration
# nodes fixed so variation across forms reflects the item-generation scheme.
s1_evaluate_cell <- function(recipe, shape, pool, pool_data = NULL,
                             n_draw = 500L, M = 20000L, seed = 1L,
                             latent_seed = NULL) {
  if (is.null(latent_seed)) {
    latent_seed <- 900000L + as.integer(strtoi(substr(
      digest::digest(shape, algo = "sha1"), 1L, 6L
    ), 16L)) %% 90000L
  }
  set.seed(seed)
  lat <- IRTsimrel::sim_latentG(M, shape = shape, seed = latent_seed)
  theta <- as.numeric(lat$theta)

  integrable <- !identical(shape, "heavy_tail")

  args <- recipe_to_args(recipe, pool, pool_data)
  out <- vector("list", n_draw)
  for (r in seq_len(n_draw)) {
    args$seed <- as.integer((as.numeric(seed) * 1000 + r) %% 2147483000)

    ip <- sim_item_params_safe(args[setdiff(names(args), "seed")], args$seed)
    d <- ip$data

    lam <- d$lambda
    g <- if (recipe$model == "3pl") d$guessing else NULL

    rr <- rho_three(1, theta, d$beta, lam,
      guessing = g, theta_var = 1,
      tail_integrable = integrable
    )
    out[[r]] <- data.frame(
      recipe_id = recipe$recipe_id, group = recipe$group, model = recipe$model,
      latent_seed = latent_seed,
      n_items = recipe$n_items, shape = shape, pool = pool, draw = r,
      rho_tilde = rr$rho_tilde, rho_psd = rr$rho_psd, w_bar = rr$w_bar,
      w_bar_truncated = rr$w_bar_truncated,
      mean_info = rr$mean_info, msem = rr$msem,
      beta_sd = stats::sd(d$beta), beta_mean = mean(d$beta),
      lambda_mean = mean(lam), lambda_sd = stats::sd(lam),
      lambda_p10 = stats::quantile(lam, 0.10, names = FALSE),
      lambda_p90 = stats::quantile(lam, 0.90, names = FALSE),
      guessing_mean = if (is.null(g)) 0 else mean(g),
      alignment = mean(d$beta) - mean(theta),
      status = rr$status,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

s1_cell_table <- function(recipes, shapes, pools) {
  cells <- expand.grid(
    recipe_id = recipes$recipe_id, shape = shapes,
    pool = pools, stringsAsFactors = FALSE
  )
  cells <- merge(cells, recipes, by = "recipe_id", sort = FALSE)
  cells <- cells[order(cells$recipe_id, cells$shape, cells$pool), ]

  cells$cond_id <- vapply(seq_len(nrow(cells)), function(i) {
    substr(
      digest::digest(paste(cells$recipe_id[i], cells$shape[i],
        cells$pool[i],
        sep = "|"
      ), algo = "sha1"),
      1L, 10L
    )
  }, character(1))
  stopifnot(!anyDuplicated(cells$cond_id))
  rownames(cells) <- NULL
  cells
}
