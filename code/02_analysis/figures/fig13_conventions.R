# Rebuild a manuscript figure from retained results.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

source("code/02_analysis/figures/fig_common.R")
inputs <- c(
  rates = cb("tables", "T_s5_form_cluster_rates.csv"),
  failures = cb("tables", "T_s5_form_cluster_mh_failures.csv"),
  bootstrap = cb("tables", "T_s5_form_cluster_bootstrap_provenance.csv")
)
expected <- c(
  rates = "3742831721b605c3ed1301037f6e526ab601bda79d88ebf73dc8a8f066ae79d7",
  failures = "1529fc57d4ad389a2d0a158cdaf1a75dc8f59f0c6804108eae17193d48410128",
  bootstrap = "5b8438ca28b24a547156e2dbe40474133f90eb1b1d7f1f6e4f3a1ad74748e768"
)
sha256 <- function(path) sub(" .*", "", system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE))
stopifnot(identical(vapply(inputs, sha256, character(1)), expected))
fcr <- read_locked(inputs["rates"]) |> filter(estimand == "successful_only")
fail <- read_locked(inputs["failures"])
boot <- read_locked(inputs["bootstrap"])
stopifnot(
  nrow(fcr) == 48L, nrow(fail) == 48L, nrow(boot) == 1L,
  boot$bootstrap_replications == 2000L,
  boot$resampling_unit == "form_cluster_block_within_structural_cell",
  boot$covariance_preserved == "all_arms_and_all_variants",
  all(fcr$n_form_clusters_planned == 3200L), all(fcr$n_structural_cells == 32L),
  all(fcr$n_form_clusters_unreachable[fcr$arm == "target_090"] == 15L),
  all(fcr$n_form_clusters_unreachable[fcr$arm != "target_090"] == 0L),
  all(is.finite(fcr$ci_low_cell_stratified_bootstrap)),
  all(is.finite(fcr$ci_high_cell_stratified_bootstrap)),
  all(fcr$ci_low_cell_stratified_bootstrap <= fcr$estimate),
  all(fcr$ci_high_cell_stratified_bootstrap >= fcr$estimate),
  max(abs(fcr$estimate - fcr$n_mh_rejected / fcr$n_mh_successful)) < 1e-12,
  sum(fcr$n_mh_attempted) == 3069600, sum(fcr$n_mh_successful) == 3057047,
  sum(fcr$n_mh_failed) == 12553
)
chk <- inner_join(fcr, fail, by = c("variant", "arm", "item_kind"), suffix = c("_rate", "_fail"))
stopifnot(
  nrow(chk) == 48L,
  all(chk$n_mh_successful_rate == chk$n_mh_successful_fail),
  all(chk$n_mh_attempted_rate == chk$n_mh_attempted_fail),
  all(chk$n_mh_failed_rate == chk$n_mh_failed_fail),
  all(chk$n_form_clusters_unreachable == chk$n_forms_unreachable),
  all(chk$n_form_clusters_contributing == 3200 - chk$n_forms_unreachable -
    chk$n_reached_forms_with_all_mh_tests_failed)
)

impact_null <- fcr |> filter(variant == "impact", item_kind == "null_clean")
stopifnot(nrow(impact_null) == 4L, all(impact_null$ci_low_cell_stratified_bootstrap > .05))
null_pair <- fcr |>
  filter(variant %in% c("baseline", "logodds"), item_kind == "null_clean") |>
  select(variant, arm, estimate, ci_low_cell_stratified_bootstrap, ci_high_cell_stratified_bootstrap) |>
  pivot_wider(names_from = variant, values_from = c(estimate, ci_low_cell_stratified_bootstrap, ci_high_cell_stratified_bootstrap))
stopifnot(
  all(null_pair$estimate_baseline == null_pair$estimate_logodds),
  all(null_pair$ci_low_cell_stratified_bootstrap_baseline == null_pair$ci_low_cell_stratified_bootstrap_logodds),
  all(null_pair$ci_high_cell_stratified_bootstrap_baseline == null_pair$ci_high_cell_stratified_bootstrap_logodds)
)
kind_map <- c(
  power = "A  Planted-item detection",
  dif_form_clean_items = "B  Clean-item flags\non DIF-containing forms",
  null_clean = "C  Fully DIF-free forms"
)
var_labs <- c(
  baseline = "Fixed difficulty shift (baseline)",
  logodds = "Fixed linear-predictor shift",
  purified = "Two-stage purification",
  impact = "Focal mean +.5 (unpurified)"
)
arm_levels <- c("uncontrolled", "target_060", "target_075", "target_090")
d <- fcr |> mutate(
  panel = factor(kind_map[item_kind], unname(kind_map)),
  var_f = factor(var_labs[variant], unname(var_labs)),
  arm_x = match(arm, arm_levels),
  plot_x = arm_x + c(-.15, -.05, .05, .15)[match(variant, names(var_labs))]
)
stopifnot(!anyNA(d$panel), !anyNA(d$var_f), !anyNA(d$plot_x))
nominal <- data.frame(panel = factor(kind_map[["null_clean"]], unname(kind_map)), y = .05)
fig13 <- ggplot(d, aes(plot_x, estimate, color = var_f, shape = var_f)) +
  geom_vline(xintercept = 1.5, linewidth = .25, color = PAL$grid) +
  geom_hline(
    data = nominal, aes(yintercept = y), inherit.aes = FALSE,
    color = PAL$ink3, linetype = "22", linewidth = .35
  ) +
  geom_line(data = d |> filter(arm != "uncontrolled"), aes(group = var_f), linewidth = .45) +
  geom_errorbar(aes(ymin = ci_low_cell_stratified_bootstrap, ymax = ci_high_cell_stratified_bootstrap),
    width = .16, linewidth = .45
  ) +
  geom_point(size = 1.5, stroke = .45, fill = "white") +
  facet_wrap(~panel, nrow = 1, scales = "free_y") +
  scale_color_manual(NULL, values = setNames(c(PAL$ink, PAL$vermillion, PAL$green, PAL$purple), unname(var_labs))) +
  scale_shape_manual(NULL, values = setNames(c(21, 22, 24, 23), unname(var_labs))) +
  scale_x_continuous(
    breaks = 1:4, labels = c("c = 1", ".60", ".75", ".90"),
    name = "Conventional arm | target index under calibration reference G",
    expand = expansion(add = .3)
  ) +
  scale_y_continuous(
    labels = label_percent(accuracy = 1), limits = c(0, NA),
    name = "MH rejection rate", expand = expansion(mult = c(0, .04))
  ) +
  labs(caption = "Points and 95% intervals: same-run pooled successful-test rates, 2,000 whole-form bootstrap samples.\nResponse abilities are normal; focal mean is +.5 only for impact. Normal/skewed calibration references are pooled.\nDashed line: .05 in fully null panel. Fifteen unavailable target-.90 arms are distinct from failed attempted tests.") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_hs() +
  theme(
    legend.position = "bottom", legend.justification = "left", legend.margin = margin(t = 1, b = 0),
    legend.key.width = unit(5, "mm"), legend.box.spacing = unit(1, "mm"),
    panel.spacing.x = unit(5, "mm"),
    strip.text = element_text(size = BASE_SIZE + .1, face = "bold", hjust = 0, margin = margin(b = 4)),
    plot.caption = element_text(size = BASE_SIZE - .8, hjust = 0, color = PAL$ink2, margin = margin(t = 3))
  )
save_fig(fig13, "fig13_conventions", MM_FULL, 101)
stopifnot(identical(vapply(inputs, sha256, character(1)), expected))
message("48 same-run rates/CIs and failure denominators verified; all three input hashes unchanged.")
print(as.data.frame(d |> select(variant, arm, item_kind, estimate, n_form_clusters_contributing)), row.names = FALSE)
