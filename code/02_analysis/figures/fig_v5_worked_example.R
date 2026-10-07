# Rebuild a manuscript figure from retained results.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

source("code/02_analysis/figures/fig_common.R")

form <- read_locked(cb("software", "normal_item_form.csv"))
calib <- read_locked(cb("software", "normal_calibration.csv"))
stopifnot(
  nrow(form) == 20L, nrow(calib) == 1L, calib$package_version == "0.3.1",
  all(abs(form$discrimination_scaled - calib$multiplier * form$discrimination_base) < 1e-9)
)
c_star <- calib$multiplier
set.seed(20260826)
nodes <- rnorm(20000)
v_hat <- var(nodes)
rho_M <- function(cc) {
  A <- mean(test_info(nodes, form$difficulty, cc * form$discrimination_base))
  v_hat * A / (1 + v_hat * A)
}
cs <- exp(seq(log(0.1), log(10), length.out = 300))
obj <- data.frame(c = cs, rho = vapply(cs, rho_M, numeric(1)))
root_here <- uniroot(function(cc) rho_M(cc) - 0.85, c(0.5, 3), tol = 1e-10)$root

message(sprintf(
  "retained c* = %.6f; recomputed root on these nodes = %.6f (v_hat = %.5f)",
  c_star, root_here, v_hat
))
stopifnot(abs(root_here - c_star) < 0.02, sum(diff(obj$rho) > 0) == length(cs) - 1L)

pA <- ggplot(obj, aes(c, rho)) +
  geom_hline(yintercept = .85, linetype = "22", color = PAL$ink3, linewidth = .35) +
  geom_line(color = PAL_FUN[["rho_tilde"]], linewidth = .7) +
  annotate("segment",
    x = c_star, xend = c_star, y = 0, yend = .85, linetype = "13",
    color = PAL$ink2, linewidth = .35
  ) +
  annotate("point", x = c_star, y = .85, color = PAL$ink, size = 1.9) +
  annot(0.105, .905, "target~widetilde(rho)*'*' == '.85'", parse = TRUE, hjust = 0, size = 2.6) +
  annot(c_star * 1.12, .06, "c^'*' == 1.278", parse = TRUE, hjust = 0, size = 2.6, color = PAL$ink) +
  scale_x_log10(
    breaks = c(0.1, 0.3, 1, 3, 10), labels = c("0.1", "0.3", "1", "3", "10"),
    name = expression("common multiplier" ~ italic(c) ~ "(log scale)")
  ) +
  scale_y_continuous(
    limits = c(0, 1), breaks = c(0, .25, .5, .75, 1), labels = lab_nolead,
    name = expression("fixed-node objective" ~ hat(widetilde(rho))[M](c))
  ) +
  labs(title = "A  The fixed-node objective") +
  theme_hs() +
  theme(plot.title = element_text(
    size = BASE_SIZE + .2, face = "bold", hjust = 0,
    color = PAL$ink, family = BASE_FAMILY, margin = margin(b = 2)
  ))

bank <- form |> transmute(beta = difficulty, base = discrimination_base, cal = discrimination_scaled)
G_df <- data.frame(theta = seq(-4, 4, .02)) |> mutate(dens = dnorm(theta))
pB <- ggplot() +
  geom_area(data = G_df, aes(theta, dens * 2.2), fill = PAL$band, color = NA) +
  geom_segment(
    data = bank, aes(x = beta, xend = beta, y = base, yend = cal), color = PAL$blue,
    linewidth = .5, arrow = arrow(length = unit(1.3, "mm"), type = "closed")
  ) +
  geom_point(data = bank, aes(beta, base), color = PAL$ink3, size = 1.3) +
  annot(3.9, 2.08, "grey: baseline slope (c = 1)", size = 2.4, hjust = 1, color = PAL$ink3) +
  annot(3.9, 1.93, "arrow: calibrated (c = 1.278)", size = 2.4, hjust = 1, color = PAL$blue) +
  scale_x_continuous(
    limits = c(-4, 4), breaks = c(-2, 0, 2),
    name = expression(theta ~ "(item difficulty" ~ beta[i] * ")")
  ) +
  scale_y_continuous(limits = c(0, 2.2), name = expression("discrimination" ~ lambda[i])) +
  labs(title = "B  Baseline and calibrated form") +
  theme_hs() +
  theme(plot.title = element_text(
    size = BASE_SIZE + .2, face = "bold", hjust = 0,
    color = PAL$ink, family = BASE_FAMILY, margin = margin(b = 2)
  ))

vals <- c(
  rho_tilde = calib$mean_information_index, rho_psd = calib$pointwise_information_index,
  w_bar = calib$harmonic_information_index
)
lad <- data.frame(fun = factor(names(vals), levels = rev(names(vals))), value = unname(vals)) |>
  mutate(lab = paste0(LAB_FUN[as.character(fun)], "*' = '*'", sub("^0\\.", ".", sprintf("%.4f", value)), "'"))
pC <- ggplot(lad, aes(value, fun, color = fun)) +
  geom_vline(xintercept = .85, linetype = "22", color = PAL$ink3, linewidth = .35) +
  geom_segment(aes(x = .80, xend = value, yend = fun), linewidth = .35, linetype = "13") +
  geom_point(size = 2.4) +
  geom_text(aes(label = lab),
    parse = TRUE, nudge_y = .32, size = 2.6, family = BASE_FAMILY,
    color = PAL$ink, hjust = .6
  ) +
  annot(.851, .55, "target", size = 2.45, hjust = 0, color = PAL$ink3) +
  scale_color_manual(values = unname(PAL_FUN[names(vals)]) |> setNames(names(vals)), guide = "none") +
  scale_x_continuous(
    limits = c(.80, .88), breaks = seq(.80, .88, .02), labels = lab_nolead,
    name = "index on 200,000 fresh draws (v = 1)"
  ) +
  scale_y_discrete(name = NULL, expand = expansion(add = c(.6, .9))) +
  labs(title = "C  Three indices at c*") +
  theme_hs() +
  theme(
    axis.text.y = element_blank(), panel.grid.major.y = element_blank(),
    plot.title = element_text(
      size = BASE_SIZE + .2, face = "bold", hjust = 0,
      color = PAL$ink, family = BASE_FAMILY, margin = margin(b = 2)
    )
  )

fig <- pA | pB | pC
save_fig(fig, "fig4_worked_example", MM_FULL, 62)
write.csv(
  data.frame(
    quantity = c("c_star_retained", "root_recomputed", "v_hat_recomputed", names(vals)),
    value = c(c_star, root_here, v_hat, unname(vals))
  ),
  file.path(DERIVED_DIR, "fig4_worked_example.csv"),
  row.names = FALSE
)
