# Rebuild a manuscript figure from retained results.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

source("code/02_analysis/figures/fig_common.R")

cost <- read_locked(ms_data("T_s1_equalization_cost.csv")) |>
  filter(reachable)
smry <- read_locked(ms_data("T_s1_equalization_summary.csv"))
stopifnot(
  nrow(cost) == 400L, all(is.finite(cost$lambda_median)),
  all(cost$lambda_median > 0)
)
check <- cost |>
  group_by(target, n_items) |>
  summarise(recomputed = median(lambda_median), n = n(), .groups = "drop") |>
  left_join(smry, by = c("target", "n_items"))
stopifnot(
  all(abs(check$recomputed - check$median_lambda_median) < 1e-12),
  all(check$n == check$n_cells)
)

BAND <- c(0.30, 2.50)
PAL_T <- c(`.70` = PAL$blue, `.85` = PAL$vermillion)

cost <- cost |> mutate(t_lab = factor(
  sprintf("%.2f", target),
  c("0.70", "0.85"), c(".70", ".85")
))
smry <- smry |> mutate(t_lab = factor(
  sprintf("%.2f", target),
  c("0.70", "0.85"), c(".70", ".85")
))

ref_anchor <- smry |> filter(t_lab == ".70", n_items == 60)
ref <- data.frame(n_items = c(5, 60)) |>
  mutate(lambda = ref_anchor$median_lambda_median *
    (n_items / 60)^(-0.5))

fig6 <- ggplot(cost, aes(n_items, lambda_median)) +
  annotate("rect",
    xmin = 4.5, xmax = 67, ymin = BAND[1], ymax = BAND[2],
    fill = PAL$band
  ) +
  geom_point(aes(color = t_lab, shape = t_lab),
    alpha = 0.25, size = 0.9,
    position = position_jitter(width = 0.012, height = 0, seed = 1)
  ) +
  geom_line(
    data = ref, aes(n_items, lambda), linetype = "31",
    color = PAL$ink3, linewidth = 0.45
  ) +
  geom_line(
    data = smry, aes(n_items, median_lambda_median, color = t_lab),
    linewidth = 0.7
  ) +
  geom_point(
    data = smry, aes(n_items, median_lambda_median, color = t_lab, shape = t_lab),
    size = 2.1
  ) +
  annot(5.2, 0.72, "reference slope -1/2",
    size = 2.45, hjust = 0, vjust = 1,
    color = PAL$ink3
  ) +
  annot(5.15, 4.6, "target~widetilde(rho)*'*' == '.85'",
    parse = TRUE,
    hjust = 0, size = 2.7, color = PAL_T[[".85"]]
  ) +
  annot(5.15, 1.28, "target~widetilde(rho)*'*' == '.70'",
    parse = TRUE,
    hjust = 0, size = 2.7, color = PAL_T[[".70"]]
  ) +
  annot(5.2, 0.34, "study diagnostic band: .30 to 2.50",
    hjust = 0, size = 2.45, color = PAL$ink3, fontface = "italic"
  ) +
  scale_color_manual(values = PAL_T, guide = "none") +
  scale_shape_manual(values = c(`.70` = 16, `.85` = 17), guide = "none") +
  scale_x_log10(
    breaks = c(5, 10, 15, 20, 40, 60),
    name = "test length (log scale)"
  ) +
  scale_y_log10(
    breaks = c(0.3, 0.5, 1, 2, 2.5, 5, 10),
    labels = c("0.3", "0.5", "1", "2", "2.5", "5", "10"),
    name = "calibrated median discrimination (log scale)"
  ) +
  theme_hs()

save_fig(fig6, "fig6_equalization", MM_MID, 78)
