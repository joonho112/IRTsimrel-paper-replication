# Rebuild a manuscript figure from retained results.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

source("code/02_analysis/figures/fig_common.R")

ms <- read_locked(cb("results", "s3", "m_sensitivity.csv")) |> filter(status == "ok")
stopifnot(nrow(ms) == 1680L, all(ms$target_rho == 0.7))
agg <- ms |>
  group_by(M) |>
  summarise(
    n = n(), holdout = mean(abs(holdout_error)),
    holdout_sd = sqrt(mean((holdout_error - mean(holdout_error))^2)),
    residual = mean(abs(calibration_residual)), .groups = "drop"
  )
stopifnot(all(agg$n == 240L), nrow(agg) == 7L)
cells <- ms |>
  group_by(M, model, latent_shape, n_items) |>
  summarise(holdout = mean(abs(holdout_error)), n = n(), .groups = "drop")
stopifnot(all(cells$n == 20L), nrow(cells) == 84L)
write.csv(agg, file.path(DERIVED_DIR, "fig5_holdout_by_M.csv"), row.names = FALSE)
write.csv(cells, file.path(DERIVED_DIR, "fig5_holdout_by_cell.csv"), row.names = FALSE)

lines <- bind_rows(
  agg |> transmute(M, err = holdout, layer = "fresh-draw error (200,000 independent draws)"),
  agg |> transmute(M, err = residual, layer = "solver residual (the calibration's own nodes)")
) |>
  mutate(layer = factor(layer, c(
    "fresh-draw error (200,000 independent draws)",
    "solver residual (the calibration's own nodes)"
  )))
labs_df <- agg |>
  filter(M %in% c(500, 20000)) |>
  mutate(lab = sub("^0\\.", ".", formatC(holdout, format = "f", digits = 4)))
Ms <- sort(unique(agg$M))

fig <- ggplot() +
  geom_point(data = cells, aes(M, holdout), color = PAL$sky, alpha = .55, size = 1.1) +
  geom_line(data = lines, aes(M, err, color = layer), linewidth = .65) +
  geom_point(data = lines, aes(M, err, color = layer, shape = layer), size = 1.9) +
  geom_text(
    data = labs_df, aes(M, holdout, label = lab), nudge_y = .28, size = 2.6,
    family = BASE_FAMILY, color = PAL$blue
  ) +
  scale_color_manual(NULL, values = c(PAL$blue, PAL$vermillion) |> setNames(levels(lines$layer))) +
  scale_shape_manual(NULL, values = c(16, 15) |> setNames(levels(lines$layer))) +
  scale_x_log10(
    breaks = Ms, labels = c("500", "1k", "2.5k", "5k", "10k", "20k", "50k"),
    name = expression("calibration node count" ~ italic(M) ~ "(log scale)")
  ) +
  scale_y_log10(
    breaks = c(1e-9, 1e-7, 1e-5, 1e-3, 1e-1),
    labels = c(
      expression(10^-9), expression(10^-7), expression(10^-5),
      expression(10^-3), expression(10^-1)
    ),
    name = "mean absolute error (log scale)"
  ) +
  coord_cartesian(ylim = c(5e-10, 2e-1)) +
  theme_hs() +
  theme(
    legend.position = c(.72, .55), legend.background = element_blank(),
    legend.key.height = unit(4, "mm")
  )
save_fig(fig, "fig5_holdout", MM_MID, 72)
print(as.data.frame(agg), digits = 4)
