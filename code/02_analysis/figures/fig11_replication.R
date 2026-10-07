# Rebuild a manuscript figure from retained results.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

source("code/02_analysis/figures/fig_common.R")
arma <- read_locked(ms_data("T_s6_arm_a.csv")) |>
  filter(te_scale == "theta_fixed", method == "Constant")
armb <- read_locked(ms_data("T_s6_arm_b_declared_primary.csv"))
stopifnot(
  nrow(arma) == 36L, all(arma$n_reps == 1000),
  all(is.finite(arma$se_ratio)), all(is.finite(arma$se_ratio_mcse)),
  nrow(armb) == 27L, all(armb$analysis == "retry_resolved_chibar"),
  all(armb$is_primary), all(armb$p_value_method == "equal_mix_chisq1_chisq2"),
  all(armb$denominator_rule == "observed"),
  sum(armb$n_observed) == 5371L, sum(armb$n_missing) == 29L,
  all(armb$n_planned == 200L),
  all(abs(armb$detection_rate - armb$n_detected / armb$n_observed) < 1e-12)
)
PAL_SD <- c(`0` = PAL$ink3, `0.2` = PAL$sky, `0.4` = PAL$blue)
SHAPE_SD <- c(`0` = 1, `0.2` = 17, `0.4` = 15)
arma <- arma |>
  mutate(
    sd_f = factor(te_sd, names(PAL_SD)),
    mean_f = factor(
      te_mean, c(0, 0.4),
      c("Reference mean 0", "Reference mean .4")
    )
  )
pa <- ggplot(arma, aes(target_rho, se_ratio,
  color = sd_f, shape = sd_f,
  group = sd_f
)) +
  geom_hline(yintercept = 1, color = PAL$ink3, linewidth = .3, linetype = "22") +
  geom_errorbar(
    aes(
      ymin = se_ratio - se_ratio_mcse,
      ymax = se_ratio + se_ratio_mcse
    ),
    width = .005, linewidth = .35, alpha = .75
  ) +
  geom_line(linewidth = .5) +
  geom_point(size = 1.6, stroke = .55) +
  facet_wrap(~mean_f, nrow = 1) +
  scale_color_manual("Reference SD", values = PAL_SD) +
  scale_shape_manual("Reference SD", values = SHAPE_SD) +
  scale_x_continuous(
    breaks = c(.5, .6, .7, .8, .9), labels = lab_nolead,
    name = "Target mean-information index"
  ) +
  scale_y_continuous(
    limits = c(.81, 1.10), breaks = c(.85, .90, .95, 1, 1.05, 1.10),
    labels = lab_nolead,
    name = "Mean model SE / empirical SD\n(constant-effect fit)"
  ) +
  labs(title = "A  Fixed length: uncertainty for the item-population effect") +
  theme_hs() +
  theme(
    plot.title = element_text(size = BASE_SIZE, face = "bold", hjust = 0),
    legend.position = "bottom", legend.margin = margin(t = 0)
  )
PAL_LEN <- c(`10` = PAL$ink3, `20` = PAL$sky, `40` = PAL$blue)
SHAPE_LEN <- c(`10` = 1, `20` = 17, `40` = 15)
armb <- armb |>
  mutate(
    len_f = factor(n_items, c(10, 20, 40)),
    sd_panel = factor(
      te_sd, c(.1, .2, .4),
      c("Reference SD .1", "Reference SD .2", "Reference SD .4")
    )
  )
pb <- ggplot(armb, aes(target_rho, detection_rate,
  color = len_f, shape = len_f,
  group = len_f
)) +
  geom_errorbar(
    aes(
      ymin = detection_rate - detection_mcse,
      ymax = detection_rate + detection_mcse
    ),
    width = .004, linewidth = .35, alpha = .75
  ) +
  geom_line(linewidth = .5) +
  geom_point(size = 1.6, stroke = .55) +
  facet_wrap(~sd_panel, nrow = 1) +
  scale_color_manual("Items", values = PAL_LEN) +
  scale_shape_manual("Items", values = SHAPE_LEN) +
  scale_x_continuous(
    breaks = c(.80, .887, .94), labels = c(".80", ".887", ".94"),
    name = "Target mean-information index"
  ) +
  scale_y_continuous(
    limits = c(0, 1), breaks = seq(0, 1, .25), labels = lab_nolead,
    name = "Heterogeneity-detection proportion\n(observed outcomes)"
  ) +
  labs(title = "B  Crossed lengths and targets: declared primary boundary test") +
  theme_hs() +
  theme(
    plot.title = element_text(size = BASE_SIZE, face = "bold", hjust = 0),
    legend.position = "bottom", legend.margin = margin(t = 0)
  )
fig11 <- pa / pb + plot_layout(heights = c(1, 1.05))
save_fig(fig11, "fig11_replication", MM_FULL, 125)
