# Rebuild a manuscript figure from retained results.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

source("code/02_analysis/figures/fig_common.R")
inputs <- c(
  s5 = cb("tables", "T_s5_metamodel.csv"),
  s6 = cb("tables", "T_s6_metamodel_collinearity.csv"),
  s5cells = cb("tables", "T_s5_rates.csv")
)
expected <- c(
  s5 = "cfd6acfe8ecebe6197fcd5715afe6a2398ac4e59f694136d4c14031412f28cc8",
  s6 = "11962b5d7b444c6fcb3f96f7edd43ea8ed5b2bbf0f6e345e3e4960116b99d2e1",
  s5cells = "9741e843332d3eb6a0beb032779662c6b4dcb4ea009671f8f223e25113ececa7"
)
sha256 <- function(path) {
  sub(" .*", "", system2(
    "shasum", c("-a", "256", shQuote(path)),
    stdout = TRUE
  ))
}
before_sha <- vapply(inputs, sha256, character(1))
stopifnot(identical(before_sha, expected))
mm5 <- read_locked(inputs[["s5"]]) |> distinct(subset, vif_rho, r2_rho_on_factors)
mm6 <- read_locked(inputs[["s6"]])
cells <- read_locked(inputs[["s5cells"]]) |> filter(kind == "power")
stopifnot(
  nrow(mm5) == 2L, nrow(mm6) == 2L,
  nrow(cells) == 128L, sum(cells$arm == "uncontrolled") == 32L,
  all(cells$n_replications == 1000L),
  identical(
    sort(unique(cells$arm)),
    sort(c("uncontrolled", "target_060", "target_075", "target_090"))
  ),
  identical(mm6$n_cells, c(9L, 27L))
)
title_style <- element_text(
  size = BASE_SIZE + .2, face = "bold", family = BASE_FAMILY,
  color = PAL$ink, hjust = 0, margin = margin(b = 4)
)
node <- function(x, y, lab, fill = "white", col = PAL$ink3, half_height = .47) {
  list(
    annotate("rect",
      xmin = x - 1.42, xmax = x + 1.42,
      ymin = y - half_height, ymax = y + half_height,
      fill = fill, color = col, linewidth = .4
    ),
    annot(x, y, lab, size = 2.7, lineheight = 1.02)
  )
}
flow <- function(x, y, xe, ye) {
  annotate(
    "segment",
    x = x, y = y, xend = xe, yend = ye, color = PAL$ink2,
    linewidth = .45, arrow = arrow(length = unit(1.4, "mm"), type = "closed")
  )
}
diagram <- ggplot() +
  node(1.5, 2.2, "Design inputs\nTest length and model\nItem-generation scheme", half_height = .68) +
  node(4.7, 2.2, "Generate and calibrate\nbank; choose multiplier c",
    fill = PAL$band, col = PAL$blue
  ) +
  node(7.9, 2.2, "Generate responses") +
  node(11.1, 2.2, "Fit models and report\nperformance") +
  flow(2.92, 2.2, 3.25, 2.2) +
  flow(6.12, 2.2, 6.45, 2.2) +
  flow(9.32, 2.2, 9.65, 2.2) +
  node(4.7, 3.9, "Calibration reference G\nand chosen index target",
    fill = PAL$band2
  ) +
  node(7.9, 3.9, "Response population\nand effect convention",
    fill = PAL$band2
  ) +
  flow(4.7, 3.43, 4.7, 2.70) +
  flow(7.9, 3.43, 7.9, 2.70) +
  node(4.7, .5, "Information curve\nand reference index") +
  flow(4.7, 1.73, 4.7, 1.00) +
  annot(9.4, .5, "At a common target, information curves\nand item parameters can still differ.",
    size = 2.55, hjust = .5, lineheight = 1.08
  ) +
  scale_x_continuous(limits = c(0, 12.7), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(-.12, 4.55), expand = expansion(mult = 0)) +
  labs(
    title = "A  From design inputs to analysis",
    caption = "Arrows show the order of calculation."
  ) +
  theme_schematic() +
  theme(
    plot.title = title_style,
    plot.caption = element_text(
      size = BASE_SIZE - .7, hjust = 0,
      color = PAL$ink2, margin = margin(t = 3)
    )
  )
coll <- bind_rows(
  mm5 |> transmute(
    row = ifelse(subset == "uncontrolled arm",
      "DIF: conventional (32 cells)", "DIF: pooled arms (128)"
    ),
    expanded = subset == "all arms pooled",
    r2 = r2_rho_on_factors, vif = vif_rho
  ),
  mm6 |> transmute(
    row = ifelse(grepl("diagonal", design),
      "Treatment: diagonal (9)", "Treatment: crossed (27)"
    ),
    expanded = !grepl("diagonal", design),
    r2 = r2_reliability_on_log_items, vif = vif
  )
) |>
  mutate(
    row = factor(row, levels = rev(c(
      "DIF: conventional (32 cells)",
      "DIF: pooled arms (128)",
      "Treatment: diagonal (9)",
      "Treatment: crossed (27)"
    ))),
    label_r2 = sprintf("%.3f", r2),
    label_vif = ifelse(vif > 10, sprintf("%.1f", vif), sprintf("%.2f", vif))
  )
stopifnot(
  nrow(coll) == 4L,
  all(abs(coll$vif - 1 / (1 - coll$r2)) < 1e-9)
)
pb <- ggplot(coll, aes(r2, row, color = expanded, shape = expanded)) +
  geom_point(size = 2.0) +
  geom_text(aes(x = r2 + .055, label = label_r2),
    hjust = 0,
    size = 2.65, family = BASE_FAMILY, color = PAL$ink
  ) +
  scale_x_continuous(
    limits = c(-.035, 1.30), breaks = c(0, .5, 1),
    labels = c("0", ".5", "1"),
    name = expression(R^2 ~ "of index on specified predictors")
  ) +
  scale_y_discrete(NULL) +
  scale_color_manual(
    values = c("FALSE" = PAL$vermillion, "TRUE" = PAL$blue),
    guide = "none"
  ) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 15), guide = "none") +
  labs(title = "B  Index predicted by design factors") +
  theme_hs() +
  theme(
    plot.title = title_style, panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = BASE_SIZE - .6)
  )
pc <- ggplot(coll, aes(vif, row, color = expanded, shape = expanded)) +
  geom_point(size = 2.0) +
  geom_text(aes(x = vif * 1.20, label = label_vif),
    hjust = 0,
    size = 2.65, family = BASE_FAMILY, color = PAL$ink
  ) +
  scale_x_log10(
    limits = c(.8, 190), breaks = c(1, 10, 100),
    name = "VIF = 1 / (1 - R-squared), log scale"
  ) +
  scale_y_discrete(NULL, labels = NULL) +
  scale_color_manual(
    values = c("FALSE" = PAL$vermillion, "TRUE" = PAL$blue),
    guide = "none"
  ) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 15), guide = "none") +
  labs(title = "C  Variance inflation") +
  theme_hs() +
  theme(plot.title = title_style, panel.grid.major.y = element_blank())
figure <- wrap_elements(full = diagram) /
  ((pb | pc) + plot_layout(widths = c(1.35, 1))) +
  plot_layout(heights = c(1.25, 1))
save_fig(figure, "fig5_confounding", MM_FULL, 114)
stopifnot(identical(vapply(inputs, sha256, character(1)), before_sha))
message("Three frozen inputs unchanged; retained diagnostics only, no regressions fitted.")
print(as.data.frame(coll), row.names = FALSE)
