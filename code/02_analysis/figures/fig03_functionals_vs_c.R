# Rebuild a manuscript figure from retained results.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

source("code/02_analysis/figures/fig_common.R")

inputs <- c(
  curves = cb("results", "theory", "curves_three_functionals.csv"),
  turns = cb("results", "theory", "turning_points.csv")
)
sha256 <- function(path) {
  sub(" .*", "", system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE))
}
before_sha <- vapply(inputs, sha256, character(1))
stopifnot(identical(unname(before_sha), c(
  "31977309ed4e20690e1a5b7229e2f06c0282bdea350a262b98eff27bf70bb36f",
  "1f4d6af2487b553f049a2c45b585b4cc29df90742f4de76a3a53fdb5686216e0"
)))
cv <- read_locked(inputs[["curves"]])
tp <- read_locked(inputs[["turns"]])

I_SHOW <- 20
cv <- cv |>
  filter(n_items == I_SHOW, guessing %in% c(0, 0.2)) |>
  mutate(
    shape_f = shape_lab(shape),
    model_f = factor(
      ifelse(guessing == 0, "2PL",
        "3PL (g = .20)"
      ),
      c("2PL", "3PL (g = .20)")
    )
  )
stopifnot(
  nrow(cv) == 728L, all(cv$n_nodes == 4001L),
  all(cv$theta_var == 1), all(cv$theta_var_source == "supplied")
)

long <- cv |>
  mutate(w_plot = ifelse(shape == "heavy_tail", w_bar_truncated, w_bar)) |>
  select(c, shape, shape_f, model_f, rho_tilde, rho_psd, w_plot) |>
  pivot_longer(c(rho_tilde, rho_psd, w_plot),
    names_to = "fun",
    values_to = "value"
  ) |>
  mutate(
    fun = recode(fun, w_plot = "w_bar"),
    fun = factor(fun, c("rho_tilde", "rho_psd", "w_bar")),
    dashed = fun == "w_bar" & shape == "heavy_tail"
  )

tp20 <- tp |>
  filter(n_items == I_SHOW, guessing %in% c(0, 0.2)) |>
  mutate(
    shape_f = shape_lab(shape),
    model_f = factor(
      ifelse(guessing == 0, "2PL",
        "3PL (g = .20)"
      ),
      c("2PL", "3PL (g = .20)")
    )
  )

turn_pts <- bind_rows(
  tp20 |> transmute(shape_f, model_f,
    fun = "rho_psd",
    c = c_turn_psd, value = psd_max
  ),
  tp20 |> transmute(shape_f, model_f,
    fun = "w_bar",
    c = c_turn_wbar, value = wbar_max
  )
) |>
  filter(is.finite(c), c <= 45) |>
  mutate(fun = factor(fun, c("rho_tilde", "rho_psd", "w_bar")))

turn_pts <- turn_pts |>
  filter(!(fun == "w_bar" & shape_f == "Heavy-tailed"))

zero_extension <- cv |>
  filter(shape == "heavy_tail") |>
  distinct(shape_f, model_f) |>
  mutate(value = 0)

fig3 <- ggplot(long, aes(c, value, color = fun, linetype = dashed)) +
  geom_hline(
    yintercept = 1, color = PAL$ink3, linewidth = 0.25,
    linetype = "22"
  ) +
  geom_line(linewidth = 0.55) +
  geom_hline(
    data = zero_extension, aes(yintercept = value),
    color = PAL_FUN[["w_bar"]], linewidth = .7,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = turn_pts, aes(c, value, color = fun), size = 1.5,
    shape = 21, fill = "white", stroke = 0.7, show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  facet_grid(model_f ~ shape_f) +
  scale_color_manual(
    values = PAL_FUN,
    breaks = c("rho_tilde", "rho_psd", "w_bar"),
    labels = expression(
      "Mean-information" ~ widetilde(rho),
      "Pointwise" ~ rho[PW], "Inverse-information" ~ bar(w)
    ),
    name = NULL,
    guide = guide_legend(nrow = 1, override.aes = list(linetype = "solid"))
  ) +
  scale_linetype_manual(
    values = c(`FALSE` = "solid", `TRUE` = "31"),
    guide = "none"
  ) +
  scale_x_log10(
    breaks = c(0.1, 1, 10),
    labels = c("0.1", "1", "10"),
    name = expression("common multiplier" ~ italic(c) ~ "(log scale)")
  ) +
  scale_y_continuous(
    limits = c(-.015, 1.04), breaks = c(0, 0.5, 1),
    labels = lab_nolead,
    name = "finite-grid index value"
  ) +
  labs(caption = paste0(
    "Heavy tails: dashed orange = finite-grid value; solid orange = population value (zero).\n",
    "Open circles mark finite-grid maxima."
  )) +
  coord_cartesian(clip = "off") +
  theme_hs(base_size = 9) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 9.5, color = PAL$ink),
    legend.key.width = unit(6, "mm"),
    legend.spacing.x = unit(3, "mm"),
    legend.margin = margin(0, 0, 2, 0),
    legend.box.margin = margin(0, 0, 2, 0),
    panel.spacing.x = unit(4, "mm"),
    strip.text.y = element_text(angle = 0, hjust = 0),
    plot.caption = element_text(
      size = 8, hjust = 0,
      color = PAL$ink2, margin = margin(t = 5)
    )
  )

save_fig(fig3, "fig3_functionals_vs_c", MM_FULL, 86)
stopifnot(identical(vapply(inputs, sha256, character(1)), before_sha))
message("Both frozen theory inputs unchanged; 8 cells x 91 multiplier values.")
