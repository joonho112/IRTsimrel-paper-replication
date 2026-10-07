# Rebuild a manuscript figure from retained results.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

source("code/02_analysis/figures/fig_common.R")

input <- cb("tables", "T_s5_form_cluster_factor_contrasts.csv")
sha256 <- function(path) {
  result <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  if (length(result) != 1L) stop("Cannot fingerprint input")
  sub(" .*", "", result)
}
expected_sha <- "a71e8f783fadc5f7d9ac73328a4dfaded30befe752df8109f6ba3e4ac3193b7d"
before_sha <- sha256(input)
stopifnot(identical(before_sha, expected_sha))

selected <- read_locked(input) |>
  filter(
    variant %in% c("baseline", "logodds"),
    contrast %in% c("length", "model"), estimand == "successful_only"
  )
stopifnot(
  nrow(selected) == 4L,
  !anyDuplicated(paste(selected$variant, selected$contrast)),
  all(selected$n_complete_paired_forms == 3200L),
  all(selected$level_high[selected$contrast == "length"] == "40"),
  all(selected$level_low[selected$contrast == "length"] == "20"),
  all(selected$level_high[selected$contrast == "model"] == "3pl_g20"),
  all(selected$level_low[selected$contrast == "model"] == "2pl"),
  all(abs(selected$delta_targeted - selected$delta_uncontrolled - selected$did) < 1e-12)
)

series <- c(
  delta_uncontrolled = "Conventional", delta_targeted = "Target .75",
  did = "Paired change"
)
display <- bind_rows(lapply(names(series), function(key) {
  data.frame(
    variant = selected$variant, contrast = selected$contrast,
    state = unname(series[[key]]),
    estimate = 100 * selected[[key]],
    lower = 100 * selected[[paste0(key, "_ci_low")]],
    upper = 100 * selected[[paste0(key, "_ci_high")]],
    stringsAsFactors = FALSE
  )
})) |>
  mutate(
    panel = factor(variant,
      levels = c("baseline", "logodds"),
      labels = c("A  Fixed latent difficulty", "B  Fixed linear predictor")
    ),
    state = factor(state, levels = unname(series)),
    y = ifelse(contrast == "length", 6, 2) - (as.integer(state) - 1),
    label = sprintf("%+.1f", estimate)
  )
stopifnot(
  nrow(display) == 12L,
  all(is.finite(display$estimate)), all(is.finite(display$lower)),
  all(is.finite(display$upper)), all(display$lower <= display$estimate),
  all(display$upper >= display$estimate)
)

headers <- expand.grid(
  panel = levels(display$panel), contrast = c("length", "model"),
  stringsAsFactors = FALSE
) |>
  mutate(
    panel = factor(panel, levels = levels(display$panel)),
    y = ifelse(contrast == "length", 6.85, 2.85),
    label = ifelse(contrast == "length", "Length: 40 - 20 items", "Model: 3PL - 2PL")
  )

colors <- c(
  Conventional = PAL$vermillion, `Target .75` = PAL$blue,
  `Paired change` = PAL$ink
)
fills <- c(
  Conventional = "white", `Target .75` = PAL$blue,
  `Paired change` = PAL$ink
)
shapes <- c(Conventional = 21, `Target .75` = 22, `Paired change` = 23)

plot <- ggplot(display, aes(x = estimate, y = y)) +
  geom_vline(xintercept = 0, color = PAL$ink3, linewidth = 0.35) +
  geom_hline(yintercept = 3.35, color = PAL$grid, linewidth = 0.35) +
  geom_errorbar(aes(xmin = lower, xmax = upper, color = state),
    orientation = "y", width = 0.15, linewidth = 0.5
  ) +
  geom_point(aes(shape = state, color = state, fill = state),
    size = 1.8, stroke = 0.35
  ) +
  geom_label(aes(x = upper + 1.25, label = label),
    hjust = 0,
    family = BASE_FAMILY, size = BASE_SIZE / 2.845, color = PAL$ink,
    fill = "white", linewidth = 0, label.padding = unit(0.35, "mm"),
    label.r = unit(0, "mm")
  ) +
  geom_text(
    data = headers, aes(x = -29.5, y = y, label = label),
    hjust = 0, fontface = "bold", family = BASE_FAMILY,
    size = BASE_SIZE / 2.845, color = PAL$ink, inherit.aes = FALSE
  ) +
  facet_wrap(~panel, nrow = 1) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = fills) +
  scale_shape_manual(values = shapes) +
  scale_x_continuous(
    limits = c(-30, 22), breaks = seq(-30, 20, 10),
    name = "Factor contrast or paired change (percentage points)",
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    limits = c(-0.5, 7.3), breaks = c(6, 5, 4, 2, 1, 0),
    labels = rep(c("Conventional", "Target .75", "Paired change"), 2),
    name = NULL, expand = expansion(mult = 0)
  ) +
  theme_hs() +
  theme(
    legend.position = "none", panel.grid.major.y = element_blank(),
    panel.spacing.x = unit(5, "mm"),
    axis.line.x = element_line(color = PAL$ink3, linewidth = 0.3),
    strip.text = element_text(
      size = BASE_SIZE + 0.2, face = "bold", hjust = 0,
      margin = margin(b = 4)
    ),
    plot.margin = margin(3, 5, 3, 2)
  )

save_fig(plot, "fig_main_dif", MM_FULL, 97)
stopifnot(identical(sha256(input), before_sha))
message("Frozen input SHA-256 unchanged: ", before_sha)
print(as.data.frame(display[, c("variant", "contrast", "state", "estimate", "lower", "upper")]),
  row.names = FALSE, digits = 6
)
