# Rebuild a manuscript figure from retained results.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

source("code/02_analysis/figures/fig_common.R")
inputs <- c(
  forms = cb("results", "s5_form_cluster", "form_records.csv.gz"),
  rates = cb("tables", "T_s5_form_cluster_rates.csv"),
  bootstrap = cb("tables", "T_s5_form_cluster_bootstrap_provenance.csv")
)
expected <- c(
  forms = "65aa2d1cabe126ed019425d2cbdf492d9c5e01f58ac1bbb477aa8a2a438357aa",
  rates = "3742831721b605c3ed1301037f6e526ab601bda79d88ebf73dc8a8f066ae79d7",
  bootstrap = "5b8438ca28b24a547156e2dbe40474133f90eb1b1d7f1f6e4f3a1ad74748e768"
)
sha256 <- function(path) {
  sub(" .*", "", system2(
    "shasum", c("-a", "256", shQuote(path)),
    stdout = TRUE
  ))
}
before_sha <- vapply(inputs, sha256, character(1))
stopifnot(identical(before_sha, expected))
forms <- readr::read_csv(inputs[["forms"]],
  show_col_types = FALSE, progress = FALSE,
  col_types = cols_only(
    structural_cell_id_full = col_character(), form_cluster_id_full = col_character(),
    replication = col_integer(), variant = col_character(), latent_shape = col_character(),
    arm = col_character(), arm_status = col_character(), item_kind = col_character(),
    n_mh_rejected = col_double(), n_mh_successful = col_double(),
    n_mh_attempted = col_double(), realized_rho = col_double(),
    impact_mean_focal = col_double(), purify = col_logical(),
    uncertainty_scope = col_character()
  )
) |>
  filter(variant == "baseline")
arms <- read_locked(inputs[["rates"]]) |>
  filter(variant == "baseline", estimand == "successful_only")
boot <- read_locked(inputs[["bootstrap"]])
stopifnot(
  nrow(forms) == 38400L, nrow(arms) == 12L,
  n_distinct(forms$form_cluster_id_full) == 3200L,
  n_distinct(forms$structural_cell_id_full) == 32L,
  all(forms$impact_mean_focal == 0), !any(forms$purify),
  all(forms$uncertainty_scope == "MH_primary_form_cluster_rerun"),
  nrow(boot) == 1L, boot$bootstrap_replications == 2000L,
  boot$resampling_unit == "form_cluster_block_within_structural_cell",
  boot$covariance_preserved == "all_arms_and_all_variants"
)
counts <- forms |> count(structural_cell_id_full, arm, item_kind)
stopifnot(nrow(counts) == 384L, all(counts$n == 100L))

cells <- forms |>
  group_by(structural_cell_id_full, latent_shape, arm, item_kind) |>
  summarise(
    rejected = sum(n_mh_rejected), successful = sum(n_mh_successful),
    rho = median(realized_rho[is.finite(realized_rho)]),
    .groups = "drop"
  ) |>
  mutate(rate = rejected / successful)
pooled_check <- forms |>
  group_by(arm, item_kind) |>
  summarise(
    rejected = sum(n_mh_rejected), successful = sum(n_mh_successful),
    attempted = sum(n_mh_attempted),
    rho = median(realized_rho[is.finite(realized_rho)]), .groups = "drop"
  ) |>
  inner_join(arms, by = c("arm", "item_kind"))
stopifnot(
  nrow(cells) == 384L, all(is.finite(cells$rate)), nrow(pooled_check) == 12L,
  all(pooled_check$rejected == pooled_check$n_mh_rejected),
  all(pooled_check$successful == pooled_check$n_mh_successful),
  all(pooled_check$attempted == pooled_check$n_mh_attempted),
  max(abs(pooled_check$rejected / pooled_check$successful -
    pooled_check$estimate)) < 1e-12,
  max(abs(pooled_check$rho - pooled_check$realized_rho_median)) < 1e-12,
  all(arms$n_form_clusters_planned == 3200L),
  all(arms$n_form_clusters_contributing[arms$arm == "target_090"] == 3185L),
  all(arms$n_form_clusters_unreachable[arms$arm == "target_090"] == 15L),
  all(arms$n_form_clusters_contributing[arms$arm != "target_090"] == 3200L),
  all(arms$n_form_clusters_unreachable[arms$arm != "target_090"] == 0L)
)
kind_map <- c(
  power = "A  Planted-item detection",
  dif_form_clean_items = "B  Clean-item flags\non DIF-containing forms",
  null_clean = "C  Fully DIF-free forms"
)
panelize <- function(d) {
  d |>
    mutate(
      panel = factor(kind_map[item_kind], levels = unname(kind_map)),
      arm_type = factor(ifelse(arm == "uncontrolled", "Conventional (c = 1)", "Targeted"),
        levels = c("Conventional (c = 1)", "Targeted")
      )
    )
}
cells <- panelize(cells) |>
  mutate(reference = factor(latent_shape,
    levels = c("normal", "skew_pos"),
    labels = c("Normal", "Skewed")
  ))
arms <- panelize(arms)
nominal <- data.frame(
  panel = factor(kind_map[["null_clean"]], levels = unname(kind_map)),
  y = .05
)
figure <- ggplot(cells, aes(rho, rate)) +
  geom_hline(
    data = nominal, aes(yintercept = y), inherit.aes = FALSE,
    color = PAL$ink3, linetype = "22", linewidth = .35
  ) +
  geom_point(aes(color = reference, shape = arm_type),
    size = .90,
    stroke = .25, alpha = .55
  ) +
  geom_errorbar(
    data = arms,
    aes(
      x = realized_rho_median, y = estimate,
      ymin = ci_low_cell_stratified_bootstrap,
      ymax = ci_high_cell_stratified_bootstrap
    ),
    inherit.aes = FALSE, width = .020, linewidth = .45, color = PAL$ink
  ) +
  geom_point(
    data = arms,
    aes(x = realized_rho_median, y = estimate, shape = arm_type),
    inherit.aes = FALSE, size = 1.6, fill = NA,
    stroke = .4, color = PAL$ink
  ) +
  facet_wrap(~panel, nrow = 1, scales = "free_y") +
  scale_color_manual("Calibration reference G",
    values = c(Normal = PAL$blue, Skewed = PAL$purple)
  ) +
  scale_shape_manual("Arm", values = c("Conventional (c = 1)" = 24, Targeted = 21)) +
  scale_x_continuous(
    breaks = c(.6, .75, .9), labels = c(".60", ".75", ".90"),
    name = "Median achieved index under the calibration reference G",
    expand = expansion(mult = .06)
  ) +
  scale_y_continuous(
    labels = label_percent(accuracy = 1), limits = c(0, NA),
    name = "MH rejection rate", expand = expansion(mult = c(0, .05))
  ) +
  labs(caption = "Small points: cell summaries. Large points and bars: pooled rates and 95% whole-form intervals.\nBoth response groups are standard normal; normal/skewed labels refer only to calibration G.") +
  guides(
    color = guide_legend(order = 1, override.aes = list(shape = 16, size = 2, alpha = 1)),
    shape = guide_legend(order = 2, override.aes = list(
      color = PAL$ink, fill = "white",
      size = 2.3, alpha = 1
    ))
  ) +
  theme_hs() +
  theme(
    panel.spacing.x = unit(5, "mm"), strip.text = element_text(
      size = BASE_SIZE + .1, face = "bold", hjust = 0, margin = margin(b = 4)
    ),
    legend.position = "bottom", legend.box = "vertical",
    legend.justification = "left", legend.margin = margin(t = 0, b = 0),
    legend.box.spacing = unit(1, "mm"), legend.spacing.y = unit(0, "mm"),
    plot.caption = element_text(
      size = BASE_SIZE - .8, hjust = 0,
      color = PAL$ink2, margin = margin(t = 3)
    )
  )
save_fig(figure, "fig9_dif_rates", MM_FULL, 96)
stopifnot(identical(vapply(inputs, sha256, character(1)), before_sha))
message("Same-run record/summary reconciliation passed for all 12 arm-kind rows; all inputs unchanged.")
print(arms[, c(
  "arm", "item_kind", "estimate", "ci_low_cell_stratified_bootstrap",
  "ci_high_cell_stratified_bootstrap", "n_mh_successful",
  "n_form_clusters_contributing"
)], row.names = FALSE)
