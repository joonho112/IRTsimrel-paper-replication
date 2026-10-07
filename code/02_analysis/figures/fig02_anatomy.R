# Rebuild a manuscript figure from retained results.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

source("code/02_analysis/figures/fig_common.R")

it <- read_locked(cb("results", "theory", "fig2_curves_item.csv"))
te <- read_locked(cb("results", "theory", "fig2_curves_test.csv"))
an <- read_locked(cb("results", "theory", "fig2_annotations.csv"))

stopifnot(
  nrow(an) == 3L, nrow(te) == 1203L, nrow(it) == 10025L,
  all(is.finite(te$tif)), all(te$tif > 0),
  all(abs(te$csem - 1 / sqrt(te$tif)) < 1e-11),
  all(an$rho_tilde >= an$rho_psd), all(an$rho_psd >= an$w_bar)
)
sums <- it |>
  group_by(key, theta) |>
  summarise(item_sum = sum(iif), .groups = "drop") |>
  left_join(te, by = c("key", "theta"))
stopifnot(all(abs(sums$item_sum - sums$tif) < 1e-11))

panel_levels <- c("i05_l1", "i10_l1", "i10_l2")
panel_labs <- c(
  i05_l1 = "A · 5 items, base scale",
  i10_l1 = "B · 10 items, base scale",
  i10_l2 = "C · 10 items, scale doubled"
)
relab <- function(d) mutate(d, panel_f = factor(key, panel_levels, panel_labs))
it <- relab(it)
te <- relab(te)
an <- relab(an)

G_df <- data.frame(theta = seq(-4, 4, 0.02)) |>
  mutate(dens = dnorm(theta))

strip_top <- theme(strip.text = element_text(
  size = BASE_SIZE + 0.2,
  face = "bold", hjust = 0
))
no_strip <- theme(strip.text = element_blank())
xsc <- scale_x_continuous(
  limits = c(-4, 4), breaks = c(-2, 0, 2),
  name = NULL
)

p_icc <- ggplot(it, aes(theta, icc, group = item)) +
  geom_line(color = PAL$ink3, linewidth = 0.32, alpha = 0.85) +
  facet_wrap(~panel_f, nrow = 1) +
  xsc +
  scale_y_continuous(
    breaks = c(0, 0.5, 1), labels = lab_nolead,
    name = "response\nprobability"
  ) +
  theme_hs() +
  strip_top

p_iif <- ggplot(it, aes(theta, iif, group = item)) +
  geom_line(color = PAL$ink3, linewidth = 0.32, alpha = 0.85) +
  facet_wrap(~panel_f, nrow = 1) +
  xsc +
  scale_y_continuous(
    limits = c(0, NA), breaks = c(0, 0.5, 1),
    labels = lab_nolead, name = "item\ninformation"
  ) +
  theme_hs() +
  no_strip

tif_max <- max(te$tif)
nz <- function(x) sub("^0\\.", ".", sprintf("%.2f", x))
an_tif <- bind_rows(
  an |> transmute(panel_f,
    y = 0.96,
    lab = sprintf("widetilde(rho)*' = '*'%s'", nz(rho_tilde))
  ),
  an |> transmute(panel_f,
    y = 0.80,
    lab = sprintf("rho[PW]*' = '*'%s'", nz(rho_psd))
  ),
  an |> transmute(panel_f,
    y = 0.64,
    lab = sprintf("bar(w)*' = '*'%s'", nz(w_bar))
  )
)

p_tif <- ggplot(te, aes(theta, tif)) +
  geom_area(
    data = G_df |> tidyr::crossing(distinct(te, panel_f)),
    aes(theta, dens * tif_max * 0.55), fill = PAL$band, color = NA
  ) +
  geom_line(color = PAL$ink, linewidth = 0.6) +
  geom_text(
    data = an_tif, aes(x = -3.85, y = tif_max * y, label = lab),
    parse = TRUE, hjust = 0, size = 2.5, family = BASE_FAMILY,
    color = PAL$ink
  ) +
  facet_wrap(~panel_f, nrow = 1) +
  xsc +
  scale_y_continuous(
    limits = c(0, tif_max * 1.04),
    name = "test\ninformation"
  ) +
  theme_hs() +
  no_strip

p_sem <- ggplot(te, aes(theta, csem)) +
  geom_line(color = PAL$vermillion, linewidth = 0.6) +
  facet_wrap(~panel_f, nrow = 1) +
  scale_x_continuous(
    limits = c(-4, 4), breaks = c(-2, 0, 2),
    name = expression(theta)
  ) +
  scale_y_continuous(
    breaks = c(0, 1, 2, 3),
    name = "information-based\nSEM"
  ) +
  coord_cartesian(ylim = c(0, 3)) +
  theme_hs() +
  no_strip

fig2 <- p_icc / p_iif / p_tif / p_sem
save_fig(fig2, "fig2_anatomy", MM_FULL, 120)
