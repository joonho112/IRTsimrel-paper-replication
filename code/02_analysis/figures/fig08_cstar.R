# Rebuild a manuscript figure from retained results.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

source("code/02_analysis/figures/fig_common.R")

input <- cb("results", "s2", "calibration.csv")
sha256 <- function(path) {
  sub(" .*", "", system2(
    "shasum", c("-a", "256", shQuote(path)),
    stdout = TRUE
  ))
}
before_sha <- sha256(input)
stopifnot(before_sha == "5b1de3f69a7f1dc5c650ecc1bb696fdebf2e1b2b62350467a4575191601e9b5c")
raw <- read_locked(input)
stopifnot(nrow(raw) == 748L, all(raw$package_version == "0.3.0"))
cal <- raw |> filter(eqc_status == "ok")
stopifnot(
  nrow(cal) == 743L, all(is.finite(cal$eqc_c_star)),
  all(cal$eqc_c_star >= .1 & cal$eqc_c_star <= 10)
)

model_labs <- c(
  rasch = "Rasch", `2pl` = "2PL",
  `3pl_g20` = "3PL, g = .20", `3pl_beta` = "3PL, Beta guessing"
)
cal <- cal |>
  mutate(
    model_f = factor(model_labs[model], unname(model_labs)),
    len_f = factor(
      n_items, c(10, 20, 40, 60),
      c("10 items", "20 items", "40 items", "60 items")
    )
  )
stopifnot(!anyNA(cal$model_f))

med <- cal |>
  group_by(model_f, len_f, target_rho) |>
  summarise(eqc_c_star = median(eqc_c_star), .groups = "drop")

PAL_LEN <- c(
  `10 items` = "#9CCEE8", `20 items` = "#56B4E9",
  `40 items` = "#0072B2", `60 items` = "#134B73"
)

fig8 <- ggplot(cal, aes(target_rho, eqc_c_star)) +
  geom_hline(yintercept = 1, color = PAL$ink3, linewidth = 0.3) +
  geom_point(
    position = position_jitter(width = 0.008, height = 0, seed = 3),
    size = 0.5, alpha = 0.3, color = PAL$ink3
  ) +
  geom_line(data = med, aes(color = len_f), linewidth = 0.6) +
  geom_point(data = med, aes(color = len_f, shape = len_f), size = 1.4) +
  facet_wrap(~model_f, nrow = 1) +
  scale_color_manual(NULL, values = PAL_LEN) +
  scale_shape_manual(NULL, values = c(
    `10 items` = 16, `20 items` = 17,
    `40 items` = 15, `60 items` = 18
  )) +
  scale_x_continuous(
    breaks = c(0.5, 0.7, 0.9),
    labels = c(".5", ".7", ".9"),
    name = expression("target" ~ widetilde(rho) * "*")
  ) +
  scale_y_log10(
    breaks = c(0.25, 0.5, 1, 2, 5, 10),
    labels = c(".25", ".5", "1", "2", "5", "10"),
    name = expression("calibrated multiplier" ~ italic(c) * "* (log scale)")
  ) +
  labs(caption = "Historical EQC 0.3.0: 743 successful / 748 retained production conditions.\nGrey points are individual conditions; colored symbols and lines summarize medians, not uncertainty.") +
  theme_hs() +
  theme(
    legend.position = "bottom", legend.margin = margin(t = 0),
    plot.caption = element_text(
      size = BASE_SIZE - .8, hjust = 0,
      color = PAL$ink2, margin = margin(t = 3)
    )
  )

save_fig(fig8, "fig8_cstar", MM_FULL, 86)
stopifnot(identical(sha256(input), before_sha))
message("Frozen calibration input unchanged; EQC only, no SAC/EQC ratios.")
print(as.data.frame(cal |> count(model_f, len_f)), row.names = FALSE)
