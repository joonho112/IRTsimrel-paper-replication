# Rebuild a manuscript figure from retained results.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

source("code/02_analysis/figures/fig_common.R")

draws <- read_locked(cb("results", "s1", "implicit_reliability.csv")) |>
  filter(is.finite(rho_tilde))

stopifnot(nrow(draws) == 100000L, all(table(draws$recipe_id) == 4000L))
cells <- read_locked(ms_data("T_s1_cells.csv"))
cov <- read_locked(ms_data("T_s1_coverage_gap.csv")) |>
  filter(corpus == "IRW public, dichotomous", functional == "rho_psd")
prov <- read_locked(ms_data("T_s1_atlas_provenance.csv"))

fam_levels <- c(
  "anchored", "baseline", "length", "discrimination",
  "spread", "guessing"
)
fam_labs <- c(
  anchored = "Named design anchors",
  baseline = "Baselines and package defaults",
  length = "Test-length sweep (2PL)",
  discrimination = "Fixed-discrimination sweep (2PL, 20 items)",
  spread = "Difficulty-spread sweep (2PL, 20 items)",
  guessing = "Guessing sweep (3PL, 20 items)"
)

meta <- cells |> distinct(recipe_id, label, group, model, n_items)

alike_ids <- meta |>
  filter(
    model == "2pl", n_items == 20, group != "discrimination",
    group != "anchored"
  ) |>
  pull(recipe_id)
stopifnot(length(alike_ids) == prov$n_alike)

sumry <- draws |>
  group_by(recipe_id) |>
  summarise(
    med = median(rho_tilde),
    q25 = quantile(rho_tilde, .25), q75 = quantile(rho_tilde, .75),
    p05 = quantile(rho_tilde, .05), p95 = quantile(rho_tilde, .95),
    .groups = "drop"
  ) |>
  left_join(meta, by = "recipe_id") |>
  mutate(
    family = factor(group, fam_levels),
    alike = recipe_id %in% alike_ids
  )

alike_meds <- sumry |>
  filter(alike) |>
  pull(med)
stopifnot(
  nrow(sumry) == 25L, nrow(cov) == 1L,
  cov$sample_analogue == "EAP empirical"
)

irw <- list(q25 = cov$emp_q25, q75 = cov$emp_q75, med = cov$emp_median)
XLIM <- c(0.32, 1)

fam_panel <- function(fam, show_x = FALSE, band_note = FALSE) {
  d <- sumry |>
    filter(family == fam) |>
    arrange(med) |>
    mutate(label_f = factor(label, levels = label))
  p <- ggplot(d, aes(y = label_f)) +
    annotate("rect",
      xmin = irw$q25, xmax = irw$q75, ymin = -Inf, ymax = Inf,
      fill = PAL$band, alpha = 0.95
    ) +
    geom_vline(
      xintercept = irw$med, color = PAL$sky, linewidth = 0.4,
      linetype = "22"
    ) +
    geom_segment(aes(x = p05, xend = p95, yend = label_f),
      linewidth = 0.35, color = PAL$ink3
    ) +
    geom_segment(aes(x = q25, xend = q75, yend = label_f),
      linewidth = 1.4, color = PAL$ink2
    ) +
    geom_point(aes(x = med, fill = alike, shape = alike),
      size = 2.1,
      color = "white", stroke = 0.45
    ) +
    scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 23), guide = "none") +
    scale_fill_manual(
      values = c(`FALSE` = PAL$ink, `TRUE` = PAL$blue),
      guide = "none"
    ) +
    scale_y_discrete(
      name = NULL,
      expand = expansion(add = c(0.5, if (band_note) 1.9 else 0.5))
    ) +
    labs(title = fam_labs[[fam]]) +
    theme_hs() +
    theme(
      plot.title = element_text(
        size = BASE_SIZE + 0.2, face = "bold",
        color = PAL$ink, hjust = 0,
        family = BASE_FAMILY,
        margin = margin(b = 1.5)
      ),
      axis.text.y = element_text(size = BASE_SIZE - 0.6),
      panel.grid.major.y = element_blank(),
      plot.margin = margin(1.5, 5, 1.5, 2)
    )
  if (band_note) {
    ymax <- nrow(d) + 1.3
    p <- p +
      annot((irw$q25 + irw$q75) / 2, ymax,
        "IRW empirical EAP: median and middle half",
        size = 2.45, color = "#2E86B5", fontface = "italic"
      )
  }
  if (show_x) {
    p + scale_x_continuous(
      limits = XLIM, breaks = seq(0.4, 1, 0.1), labels = lab_nolead,
      expand = expansion(mult = c(0, 0.01)),
      name = expression("Mean-information index" ~ widetilde(rho) ~ "at" ~ italic(c) == 1)
    )
  } else {
    p + scale_x_continuous(
      limits = XLIM, breaks = seq(0.4, 1, 0.1),
      labels = NULL, expand = expansion(mult = c(0, 0.01)),
      name = NULL
    )
  }
}

panels <- list(
  fam_panel("anchored", band_note = TRUE),
  fam_panel("baseline"),
  fam_panel("length"),
  fam_panel("discrimination"),
  fam_panel("spread"),
  fam_panel("guessing", show_x = TRUE)
)

hts <- c(3 + 2.6, 6 + 1, 5 + 1, 4 + 1, 4 + 1, 3 + 2.1)
fig4 <- wrap_plots(panels, ncol = 1, heights = hts)

save_fig(fig4, "fig4_atlas", MM_FULL, 132)
