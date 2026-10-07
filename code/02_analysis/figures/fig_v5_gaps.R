# Rebuild a manuscript figure from retained results.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

source("code/02_analysis/figures/fig_common.R")

cal <- read_locked(cb("results", "s2", "calibration.csv")) |> filter(eqc_status == "ok")
stopifnot(nrow(cal) == 743L)
model_labs <- c(rasch = "Rasch", `2pl` = "2PL", `3pl_g20` = "3PL, g = .20", `3pl_beta` = "3PL, Beta guessing")
d <- cal |>
  transmute(
    model_f = factor(model_labs[model], unname(model_labs)),
    shape_f = shape_lab(latent_shape), target_rho,
    gap_pw = rho_tilde_at_cstar - rho_psd_at_cstar,
    gap_w = ifelse(is.na(w_bar_at_cstar) | w_bar_at_cstar == "NA", NA_real_,
      rho_tilde_at_cstar - as.numeric(w_bar_at_cstar)
    )
  )
stopifnot(sum(!is.na(d$gap_w)) == 559L, all(d$gap_pw >= -1e-12))
long <- bind_rows(
  d |> transmute(model_f, shape_f, target_rho, gap = gap_pw, which = "pointwise"),
  d |> filter(!is.na(gap_w)) |> transmute(model_f, shape_f, target_rho, gap = gap_w, which = "inverse-information")
) |>
  mutate(which = factor(which, c("pointwise", "inverse-information")))
summ <- long |>
  group_by(which) |>
  summarise(n = n(), median = median(gap), p90 = quantile(gap, .9), max = max(gap), .groups = "drop")
print(as.data.frame(summ), digits = 4)
write.csv(summ, file.path(DERIVED_DIR, "figA1_gaps_summary.csv"), row.names = FALSE)

med <- long |>
  group_by(model_f, shape_f, which) |>
  summarise(n = n(), median = median(gap), .groups = "drop")
stopifnot(nrow(med) == 28L, all(med$n > 0L))
write.csv(med, file.path(DERIVED_DIR, "figA1_gaps_medians_by_group.csv"), row.names = FALSE)

fig <- ggplot(long, aes(shape_f, gap, color = which)) +
  geom_hline(yintercept = 0, color = PAL$ink3, linewidth = .3) +
  geom_point(
    position = position_jitterdodge(jitter.width = .25, dodge.width = .6, seed = 5),
    size = .6, alpha = .45
  ) +
  geom_point(
    data = med, aes(shape_f, median, group = which), inherit.aes = FALSE,
    position = position_dodge(width = .6), shape = 95, size = 6, colour = PAL$ink
  ) +
  facet_wrap(~model_f, nrow = 1) +
  scale_color_manual(NULL,
    values = c(
      pointwise = PAL_FUN[["rho_psd"]],
      `inverse-information` = PAL_FUN[["w_bar"]]
    ),
    labels = c(expression(widetilde(rho) - rho[PW]), expression(widetilde(rho) - bar(w)))
  ) +
  scale_y_continuous(
    name = expression("gap below" ~ widetilde(rho) ~ "at" ~ italic(c) * "* (square-root scale)"),
    trans = "sqrt", breaks = c(0, .01, .05, .1, .25, .5, .9), labels = lab_nolead
  ) +
  scale_x_discrete(name = NULL, labels = c("normal", "skewed", "bimodal", "heavy")) +
  theme_hs() +
  theme(
    legend.position = "bottom", legend.margin = margin(t = 0),
    axis.text.x = element_text(size = BASE_SIZE - 1.2, angle = 30, hjust = 1)
  )
save_fig(fig, "figA1_gaps", MM_FULL, 70)
