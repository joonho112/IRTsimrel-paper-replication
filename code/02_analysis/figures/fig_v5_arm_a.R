# Rebuild a manuscript figure from retained results.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

source("code/02_analysis/figures/fig_common.R")

arma <- read_locked(cb("tables", "T_s6_arm_a.csv")) |>
  filter(
    te_scale == "theta_fixed", method == "Constant",
    (te_mean == .4 & te_sd == 0) | (te_mean == 0 & te_sd == .4)
  ) |>
  arrange(desc(te_mean), target_rho)
targets <- c(.50, .60, .70, .80, .887, .90)
target_labels <- c(".50", ".60", ".70", ".80", ".887", ".90")
stopifnot(
  nrow(arma) == 12L, all(arma$n_reps == 1000L),
  all(arma$rung_seed_relation == "independent"),
  identical(arma$target_rho[arma$te_mean == .4], targets),
  identical(arma$target_rho[arma$te_mean == 0], targets)
)
arma$successes <- round(arma$reject_rate * arma$n_reps)
stopifnot(all(abs(arma$reject_rate * arma$n_reps - arma$successes) < 1e-9))
cis <- t(vapply(
  seq_len(nrow(arma)), function(i) {
    as.numeric(binom.test(arma$successes[i], arma$n_reps[i], conf.level = .95)$conf.int)
  },
  numeric(2)
))
arma$lower <- cis[, 1]
arma$upper <- cis[, 2]
stopifnot(sum(arma$lower[arma$te_mean == 0] > .05) == 5L)
write.csv(arma |> select(te_mean, te_sd, target_rho, successes, n_reps, reject_rate, lower, upper),
  file.path(DERIVED_DIR, "fig7_arm_a.csv"),
  row.names = FALSE
)

layer_names <- c(
  "Power increases with the information target",
  "Excess rejection persists across information targets"
)
arma <- arma |>
  mutate(
    layer = factor(ifelse(te_mean == .4, layer_names[1], layer_names[2]), levels = layer_names),
    target = factor(target_rho, levels = targets, labels = target_labels),
    label = sprintf("%.1f", reject_rate * 100),
    text_y = upper + ifelse(te_mean == .4, .030, .010)
  )
bounds <- data.frame(
  layer = factor(rep(layer_names, each = 2), levels = layer_names),
  target = factor(".50", levels = target_labels),
  y = c(.38, .84, 0, .135)
)
null_line <- data.frame(layer = factor(layer_names[2], levels = layer_names), nominal = .05)

fig <- ggplot(arma, aes(target, reject_rate, color = layer)) +
  geom_blank(data = bounds, aes(target, y), inherit.aes = FALSE) +
  geom_hline(
    data = null_line, aes(yintercept = nominal), color = PAL$ink3,
    linetype = "22", linewidth = .35, inherit.aes = FALSE
  ) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .16, linewidth = .45) +
  geom_point(size = 1.7) +
  geom_text(aes(y = text_y, label = label),
    color = PAL$ink, family = BASE_FAMILY,
    size = (BASE_SIZE - .4) / 2.845
  ) +
  facet_wrap(~layer, ncol = 1, scales = "free_y") +
  scale_color_manual(values = setNames(c(PAL$blue, PAL$vermillion), layer_names), guide = "none") +
  scale_y_continuous(
    labels = label_percent(accuracy = 1),
    breaks = function(lim) if (max(lim) > .2) seq(.4, .8, .1) else c(0, .05, .10),
    expand = expansion(mult = 0)
  ) +
  labs(
    subtitle = "Fixed across targets: 500 persons, 20 items, ability-scale effects",
    x = expression("target mean-information index" ~ widetilde(rho) * "*"),
    y = "constant-effect model rejection"
  ) +
  theme_hs() +
  theme(
    panel.grid.major.x = element_blank(), panel.spacing.y = unit(5, "mm"),
    plot.subtitle = element_text(
      size = BASE_SIZE, color = PAL$ink2, hjust = 0,
      margin = margin(b = 6)
    ),
    strip.text = element_text(size = BASE_SIZE, face = "bold", hjust = 0)
  )
save_fig(fig, "fig7_arm_a", MM_MID, 96)
print(as.data.frame(arma[, c("te_mean", "te_sd", "target_rho", "reject_rate", "lower", "upper")]), digits = 4)
