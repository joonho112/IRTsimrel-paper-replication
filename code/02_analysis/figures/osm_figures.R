# Rebuild a manuscript figure from retained results.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

source("code/02_analysis/figures/fig_common.R")

choice <- sub("^--figure=", "", grep("^--figure=", commandArgs(TRUE), value = TRUE))
stopifnot(length(choice) <= 1L, !length(choice) || choice %in% c("B1", "D1", "E1"))
selected <- function(x) !length(choice) || identical(choice, x)

if (selected("B1")) {
  ms <- read_locked(cb("results", "s3", "m_sensitivity.csv")) |>
    filter(status == "ok")
  model_labs <- c(
    `2pl` = "2PL", rasch = "Rasch", `3pl_g20` = "3PL, g = .20",
    `3pl_beta` = "3PL, Beta guessing"
  )
  msum <- ms |>
    group_by(model, M) |>
    summarise(
      residual = mean(abs(calibration_residual)),
      holdout = mean(abs(holdout_error)), .groups = "drop"
    ) |>
    pivot_longer(c(residual, holdout),
      names_to = "metric",
      values_to = "err"
    ) |>
    mutate(
      model_f = factor(
        model_labs[model],
        c(
          "Rasch", "2PL", "3PL, g = .20",
          "3PL, Beta guessing"
        )
      ),
      metric_f = factor(
        metric, c("holdout", "residual"),
        c(
          "holdout error (independent quadrature)",
          "solver residual (same quadrature)"
        )
      )
    )
  stopifnot(!anyNA(msum$model_f))

  figB1 <- ggplot(msum, aes(M, err, color = metric_f)) +
    geom_hline(
      yintercept = 5e-3, linetype = "22", color = PAL$ink3,
      linewidth = 0.3
    ) +
    annot(1.8e4, 6.6e-3, "5 %*% 10^-3",
      parse = TRUE, size = 2.4,
      color = PAL$ink3
    ) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.6) +
    facet_wrap(~model_f, nrow = 1) +
    scale_color_manual(NULL, values = c(PAL$blue, PAL$vermillion) |>
      setNames(levels(msum$metric_f))) +
    scale_x_log10(
      breaks = c(1e3, 1e4, 5e4),
      labels = c("1k", "10k", "50k"),
      name = expression("integration draws" ~ italic(M) ~ "(log scale)")
    ) +
    scale_y_log10(
      breaks = c(1e-9, 1e-6, 1e-3),
      labels = c(
        expression(10^-9), expression(10^-6),
        expression(10^-3)
      ),
      name = "mean absolute error\n(log scale)"
    ) +
    theme_hs() +
    theme(legend.position = "bottom", legend.margin = margin(t = 0))
  save_fig(figB1, "figB1_m_sensitivity", MM_FULL, 58)
}

if (selected("D1")) {
  vd <- read_locked(ms_data("T_s1_variance_decomposition.csv")) |>
    filter(factor != "Residuals")
  lab <- c(
    f_length = "Test length",
    f_lambda = "Discrimination specification",
    f_spread = "Difficulty-spread specification",
    f_guess = "Guessing specification",
    f_shape = "Latent shape",
    f_pool = "Item source",
    beta_sd = "Realized difficulty spread"
  )
  vd$factor[vd$factor == "factor(pool)"] <- "f_pool"
  vd <- vd |>
    mutate(
      factor_f = lab[factor],
      arm_f = factor(
        model_arm,
        c("parametric only", "both arms, realized spread"),
        c(
          "Parametric arm",
          "Both arms, with realized difficulty spread"
        )
      )
    )
  stopifnot(!anyNA(vd$factor_f))
  ord <- vd |>
    group_by(factor_f) |>
    summarise(m = max(share_of_explained), .groups = "drop") |>
    arrange(m)
  vd$factor_f <- factor(vd$factor_f, ord$factor_f)

  figD1 <- ggplot(vd, aes(share_of_explained, factor_f)) +
    geom_col(fill = PAL$blue, width = 0.62) +
    geom_text(aes(label = sprintf("%.2f%%", 100 * share_of_explained)),
      hjust = -0.15, size = 2.5, family = BASE_FAMILY,
      color = PAL$ink2
    ) +
    facet_wrap(~arm_f, nrow = 1) +
    scale_x_continuous(
      limits = c(0, 0.72), breaks = seq(0, 0.6, 0.2),
      labels = scales::percent_format(accuracy = 1),
      name = "share of fitted-model explained sum of squares"
    ) +
    scale_y_discrete(NULL) +
    theme_hs()
  save_fig(figD1, "figD1_decomposition", MM_FULL, 58)
}

if (selected("E1")) {
  rp <- read_locked(cb("results", "s2", "replication.csv"))
  rp <- rp |> mutate(
    shape_f = shape_lab(latent_shape),
    bias = tilde_mean - population_tilde
  )
  stopifnot(!anyNA(rp$shape_f))

  sd_med <- rp |>
    group_by(n_persons) |>
    summarise(tilde_sd = median(tilde_sd), .groups = "drop")
  ref <- data.frame(n_persons = c(100, 2000)) |>
    mutate(tilde_sd = sd_med$tilde_sd[sd_med$n_persons == 100] *
      sqrt(100 / n_persons))

  pE1 <- ggplot(rp, aes(n_persons, tilde_sd)) +
    geom_line(aes(group = cond_id),
      color = PAL$grid, linewidth = 0.2,
      alpha = 0.5
    ) +
    geom_line(
      data = ref, linetype = "31", color = PAL$vermillion,
      linewidth = 0.5
    ) +
    geom_line(data = sd_med, color = PAL$blue, linewidth = 0.8) +
    geom_point(data = sd_med, color = PAL$blue, size = 1.9) +
    annot(700, 0.045, "square-root reference",
      size = 2.5,
      color = PAL$vermillion
    ) +
    scale_x_log10(
      breaks = c(100, 250, 500, 1000, 2000),
      labels = c("100", "250", "500", "1,000", "2,000"),
      name = expression("persons drawn" ~ italic(N) ~ "(log scale)")
    ) +
    scale_y_log10(name = expression("SD of realized" ~ widetilde(rho) ~ "(log scale)")) +
    labs(title = "A · Spread across person redraws") +
    theme_hs() +
    theme(plot.title = element_text(
      size = BASE_SIZE + 0.2, face = "bold",
      color = PAL$ink, hjust = 0,
      family = BASE_FAMILY, margin = margin(b = 2)
    ))

  bias_med <- rp |>
    group_by(n_persons, shape_f) |>
    summarise(bias = median(bias), .groups = "drop")

  pE2 <- ggplot(bias_med, aes(n_persons, bias, color = shape_f)) +
    geom_hline(yintercept = 0, color = PAL$ink3, linewidth = 0.3) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.7) +
    scale_color_manual(NULL, values = PAL_SHAPE) +
    scale_x_log10(
      breaks = c(100, 250, 500, 1000, 2000),
      labels = c("100", "250", "500", "1,000", "2,000"),
      name = expression("persons drawn" ~ italic(N) ~ "(log scale)")
    ) +
    scale_y_continuous(
      name = expression("median bias of realized" ~ widetilde(rho)),
      labels = function(x) sprintf("%.3f", x)
    ) +
    labs(title = "B · Deviation from calibration-rule value") +
    theme_hs() +
    theme(
      plot.title = element_text(
        size = BASE_SIZE + 0.2, face = "bold",
        color = PAL$ink, hjust = 0,
        family = BASE_FAMILY, margin = margin(b = 2)
      ),
      legend.position = "bottom", legend.margin = margin(t = 0)
    )

  figE1 <- pE1 + pE2 + plot_layout(widths = c(1, 1))
  save_fig(figE1, "figE1_replication", MM_FULL, 66)
}
