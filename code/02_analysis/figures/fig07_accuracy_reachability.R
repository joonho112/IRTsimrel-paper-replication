# Rebuild a manuscript figure from retained results.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

source("code/02_analysis/figures/fig_common.R")

inputs <- c(
  calibration = cb("results", "s2", "calibration.csv"),
  accuracy = cb("tables", "T_s2_accuracy.csv"),
  reachability = cb("tables", "T_s2_reachability_unpruned.csv")
)
expected_sha <- c(
  calibration = "5b1de3f69a7f1dc5c650ecc1bb696fdebf2e1b2b62350467a4575191601e9b5c",
  accuracy = "ce97c041cc610866597bff8c1b37e8cf304836efb36442ba59c988b639f79a7e",
  reachability = "0173f24aa66a397038fe7c1d8ac3a849ea7ba30f54914ee551745ff1089d565d"
)
sha256 <- function(path) {
  result <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  if (length(result) != 1L) stop("Cannot fingerprint input")
  sub(" .*", "", result)
}
before_sha <- vapply(inputs, sha256, character(1))
stopifnot(identical(before_sha, expected_sha))
cal <- read_locked(inputs[["calibration"]])
accuracy <- read_locked(inputs[["accuracy"]]) |>
  filter(grouping == "algorithm", algorithm %in% c("EQC", "SAC (rho-tilde)"))
reach <- read_locked(inputs[["reachability"]])
stopifnot(
  nrow(cal) == 748L, nrow(accuracy) == 2L, nrow(reach) == 6L,
  all(cal$package_version == "0.3.0"),
  all(reach$endpoint_rule == "strict_interior_lo_lt_target_lt_hi"),
  all(reach$feasibility_quadrature_M == 20000L),
  sum(reach$tilde_eligible) == 768L,
  sum(reach$tilde_unreachable) == 20L,
  sum(reach$wbar_eligible) == 576L,
  sum(reach$wbar_unreachable) == 167L,
  sum(reach$tilde_eligible - reach$wbar_eligible) == 192L
)

acc <- bind_rows(
  cal |> filter(eqc_status == "ok") |>
    transmute(target = target_rho, delta = eqc_delta, arm = "EQC"),
  cal |> filter(sac_info_status == "ok") |>
    transmute(target = target_rho, delta = sac_info_delta, arm = "SAC")
)
stopifnot(
  sum(acc$arm == "EQC") == 743L, sum(acc$arm == "SAC") == 733L,
  all(is.finite(acc$delta)),
  accuracy$n_ok[accuracy$algorithm == "EQC"] == 743L,
  accuracy$n_ok[accuracy$algorithm == "SAC (rho-tilde)"] == 733L
)

for (arm in c("EQC", "SAC")) {
  key <- if (arm == "EQC") "EQC" else "SAC (rho-tilde)"
  vals <- abs(acc$delta[acc$arm == arm])
  row <- accuracy[accuracy$algorithm == key, ]
  stopifnot(
    abs(mean(vals) - row$mae) <= max(1e-12, abs(row$mae) * 0.005),
    abs(max(vals) - row$max_abs) <= max(1e-12, abs(row$max_abs) * 0.005)
  )
}
targets <- c(.50, .60, .70, .80, .90, .95)
target_labels <- c(".50", ".60", ".70", ".80", ".90", ".95")
acc <- acc |>
  mutate(
    target = factor(target, levels = targets, labels = target_labels),
    layer = factor(arm,
      levels = c("EQC", "SAC"),
      labels = c(
        "EQC: own objective (743/748)",
        "SAC: fresh evaluation (733/748)"
      )
    )
  )
residual_labels <- function(x) {
  vapply(x, function(v) {
    if (is.na(v)) {
      return("")
    }
    if (v == 0) {
      return("0")
    }
    if (abs(v) < 1e-6) {
      return(sub("e-0", "e-", formatC(v, format = "e", digits = 0)))
    }
    sub("^(-?)0\\.", "\\1.", sprintf("%.02f", v))
  }, character(1))
}
panel_title <- element_text(
  size = BASE_SIZE + .2, face = "bold",
  family = BASE_FAMILY, color = PAL$ink,
  hjust = 0, margin = margin(b = 5)
)
pa <- ggplot(acc, aes(target, delta)) +
  geom_hline(yintercept = 0, color = PAL$ink3, linewidth = .3) +
  geom_point(
    position = position_jitter(width = .20, height = 0, seed = 7),
    size = .55, alpha = .42, color = PAL$ink
  ) +
  facet_wrap(~layer, ncol = 1, scales = "free_y") +
  scale_y_continuous(
    labels = residual_labels,
    breaks = function(lim) {
      if (diff(lim) < 1e-6) {
        seq(-4e-8, 4e-8, 2e-8)
      } else {
        seq(-.04, .04, .02)
      }
    },
    expand = expansion(mult = .12)
  ) +
  labs(
    title = "A  Successful-cell numerical checks",
    x = "Target mean-information index", y = "Achieved - target",
    caption = "Different evaluation rules and y-axis scales"
  ) +
  theme_hs() +
  theme(
    plot.title = panel_title, plot.caption = element_text(
      size = BASE_SIZE - .8, hjust = 0, color = PAL$ink2, margin = margin(t = 5)
    ),
    panel.grid.major.x = element_blank(), panel.spacing.y = unit(5, "mm"),
    strip.text = element_text(size = BASE_SIZE, face = "bold", hjust = 0),
    plot.margin = margin(3, 5, 3, 1)
  )

reach_long <- bind_rows(
  reach |> transmute(target,
    state = "Mean info (20/768)",
    count = tilde_unreachable, eligible = tilde_eligible,
    share = tilde_share
  ),
  reach |> transmute(target,
    state = "Harmonic (167/576)",
    count = wbar_unreachable, eligible = wbar_eligible,
    share = wbar_share
  )
) |>
  mutate(
    state = factor(state, levels = c("Mean info (20/768)", "Harmonic (167/576)")),
    y = match(target, targets) + ifelse(as.integer(state) == 1, .17, -.17),
    label = paste0(count, "/", eligible)
  )
stopifnot(
  nrow(reach_long) == 12L,
  all(abs(reach_long$share - reach_long$count / reach_long$eligible) < 1e-12)
)
fill <- c(
  "Mean info (20/768)" = PAL_FUN[["rho_tilde"]],
  "Harmonic (167/576)" = PAL_FUN[["w_bar"]]
)
pb <- ggplot(reach_long, aes(share, y, fill = state)) +
  geom_col(orientation = "y", width = .25) +
  geom_text(aes(x = share + .023, label = label),
    hjust = 0,
    family = BASE_FAMILY, size = (BASE_SIZE - .3) / 2.845, color = PAL$ink
  ) +
  scale_fill_manual(values = fill, name = NULL) +
  scale_x_continuous(
    limits = c(0, 1.03), breaks = seq(0, 1, .25),
    labels = label_percent(accuracy = 1),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    limits = c(.5, 6.5), breaks = seq_along(targets),
    labels = target_labels, expand = expansion(mult = 0)
  ) +
  labs(
    title = "B  Bounded reachability census",
    x = "Unreachable among eligible cells", y = "Target",
    caption = "Harmonic: 192 polynomial-tail cells excluded\nTheir extended harmonic index is 0"
  ) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_hs() +
  theme(
    plot.title = panel_title, legend.position = "top",
    legend.justification = "left", legend.margin = margin(0, 0, 3, 0),
    legend.box.spacing = unit(0, "mm"),
    panel.grid.major.y = element_blank(),
    plot.caption = element_text(
      size = BASE_SIZE - .8, hjust = 0,
      color = PAL$ink2, margin = margin(t = 5)
    ),
    plot.margin = margin(3, 3, 3, 2)
  )
figure <- pa + pb + plot_layout(widths = c(1, 1))
save_fig(figure, "fig7_accuracy_reachability", MM_FULL, 112)
stopifnot(identical(vapply(inputs, sha256, character(1)), before_sha))
message("All three frozen-input SHA-256 values unchanged.")
print(as.data.frame(reach_long[, c("target", "state", "count", "eligible", "share")]),
  row.names = FALSE
)
