# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

# Supplementary examples share the same checkpoint runner as validation.
# Smoke mode retains the original sample sizes and numerical integration,
# but uses two response replications (or two item forms).

run_extension <- function(study, mode, pool, workers) {
  require_packages(c("digest", "jsonlite", "yaml"))
  use_irtsimrel("0.3.0")
  load_scientific_helpers()
  source(repo_path("code", "lib", "study_workers.R"))
  source(repo_path("code", "lib", "run_jobs.R"))
  source(repo_path("code", "01_simulation", "validation.R"))
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  context <- list()
  if (study == "recipes") {
    recipes <- read.csv(repo_path("config", "design_recipes.csv"))
    design <- s1_cell_table(
      recipes, c("normal", "skew_pos", "bimodal", "heavy_tail"),
      c("parametric", "irw")
    )
    design <- select_pool(design, pool)
    if (mode == "smoke") design <- design[design$recipe_id == recipes$recipe_id[1], ]
    context <- list(
      recipes = recipes, pool_data = read_irw_pool(any(design$pool == "irw")),
      n_draw = if (mode == "smoke") 2L else 500L
    )
    worker <- function(row, ctx) {
      seed <- as.integer((20260826L %% 90000L +
        strtoi(substr(row$cond_id, 1, 6), 16L)) %% 900000L + 1000L)
      s1_evaluate_cell(ctx$recipes[ctx$recipes$recipe_id == row$recipe_id, ],
        row$shape, row$pool, ctx$pool_data,
        n_draw = ctx$n_draw, M = 20000L, seed = seed
      )
    }
  } else if (study %in% c("treatment", "treatment-b")) {
    require_packages("lme4")
    cfg <- yaml::read_yaml(repo_path("config", "s6_ilhte_design.yml"))
    seed0 <- as.integer(cfg$global_seed %% 900000L)
    set.seed(seed0)
    quad <- rnorm(20000L)
    if (study == "treatment") {
      a <- cfg$arm_a
      design <- expand.grid(
        target_rho = unlist(a$reliability_ladder),
        te_mean = unlist(a$te_mean), te_sd = unlist(a$te_sd),
        te_scale = c(cfg$te_scale$primary, cfg$te_scale$sensitivity),
        rep = seq_len(a$reps), stringsAsFactors = FALSE
      )
      design <- design[design$te_scale == cfg$te_scale$primary |
        design$rep <= ceiling(a$reps / 10), ]
      design$n_sub <- a$n_sub
      design$n_items <- a$n_items
      design$cond_id <- sprintf(
        "a_%s_r%03.0f_m%02.0f_s%02.0f_%04d",
        substr(design$te_scale, 1, 5), design$target_rho * 1000,
        design$te_mean * 100, design$te_sd * 100, design$rep
      )
      if (mode == "smoke") {
        design <- design[design$rep <= 2 &
          design$te_scale == "theta_fixed" & design$target_rho %in% c(.5, .9) &
          design$te_mean == .4 & design$te_sd == .4, ]
      }
    } else {
      a <- cfg$arm_b
      design <- expand.grid(
        n_items = unlist(a$n_items),
        target_rho = unlist(a$reliability), te_sd = unlist(a$te_sd),
        rep = seq_len(a$reps), stringsAsFactors = FALSE
      )
      design$n_sub <- a$n_sub
      design$te_mean <- a$te_mean
      design$te_scale <- cfg$te_scale$primary
      design$cond_id <- sprintf(
        "b_I%02d_r%03.0f_s%02.0f_%04d", design$n_items,
        design$target_rho * 1000, design$te_sd * 100, design$rep
      )
      if (mode == "smoke") {
        design <- design[design$rep <= 2 &
          design$n_items == 10 & design$target_rho == .8 & design$te_sd == .4, ]
      }
      message("Arm B returns initial fits. See docs/reproduction-scope.md for the retained retry analysis.")
    }
    context <- list(quad = quad, seed0 = seed0, disc = cfg$original$disc)
    worker <- function(row, ctx) {
      value <- ilhte_replicate_targeted(as.integer(row$n_sub), as.integer(row$n_items),
        row$target_rho, row$te_mean, row$te_sd, row$te_scale,
        theta_quad = ctx$quad,
        seed = ctx$seed0 + as.integer(row$rep) * 19L +
          as.integer(row$target_rho * 1000) + as.integer(row$n_items) * 3L,
        disc_reference = ctx$disc
      )
      value$cond_id <- row$cond_id
      value$rep <- row$rep
      value
    }
  } else if (study == "dif") {
    base_cfg <- yaml::read_yaml(repo_path("config", "s5_dif_design.yml"))
    cfg <- yaml::read_yaml(repo_path("config", "s5_form_cluster_design.yml"))
    structural <- s5_structural_cells(base_cfg)
    # Build all named streams before selecting cells; selection cannot change a seed.
    registry <- s5_build_seed_registry(structural, cfg)
    design <- select_pool(structural, pool)
    if (mode == "smoke") {
      design <- design[design$model == "2pl" &
        design$latent_shape == "normal" & design$n_items == 20 & design$n_per_group == 500, ]
    }
    context <- list(
      base_cfg = base_cfg, cfg = cfg, rng_kind = RNGkind(),
      catalog = s5_condition_catalog(base_cfg, cfg, structural), seed_registry = registry,
      pool_data = read_irw_pool(any(design$pool == "irw")),
      M_CALIB = as.integer(cfg$quadrature$calibration_M),
      M_EVAL = as.integer(cfg$quadrature$evaluation_M),
      ROOT_INTERVAL = as.numeric(unlist(cfg$quadrature$root_interval)),
      ROOT_TOL = as.numeric(cfg$quadrature$root_tolerance),
      R_REP = if (mode == "smoke") 2L else as.integer(cfg$replications_per_structural_cell),
      VARIANTS = s5_variant_table(cfg), ARMS = as.character(unlist(cfg$arms))
    )
    worker <- original_worker("dif_worker_original")
  } else if (study == "scores") {
    require_packages("TAM")
    design <- expand.grid(delta = c(0, .2), condition = c("balanced", "imbalanced"))
    design$focal_sd <- ifelse(design$condition == "imbalanced", .6, 1)
    design$cond_id <- sprintf("s5b_d%02.0f_%s", design$delta * 100, design$condition)
    design$seed <- 910000L + seq_len(nrow(design)) * 37L
    context <- list(
      N_PER_GROUP = 500L, N_ITEMS = 20L,
      N_REPS = if (mode == "smoke") 2L else 200L, TARGET_REF = .80, ALPHA = .05
    )
    worker <- original_worker("scores_worker_original")
  } else {
    stop("Unknown supplementary study: ", study)
  }
  run_jobs(design, worker, context, study, mode, workers)
}
