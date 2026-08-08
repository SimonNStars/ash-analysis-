# Mt. Hood Ash Dispersal -- Sensitivity Analysis (R companion)
#
# Reads the OAT sensitivity CSVs written by tephra2_sensitivity.ipynb
# (oat_plume_height.csv, oat_eruption_mass.csv, oat_diffusion_coef.csv --
# columns: <parameter>, Rhododendron, Parkdale, `Govt. Camp`) and reproduces
# the same plots/summary in R, for comparison against the notebook's PNGs.
# Does not run Tephra2 itself -- purely a downstream visualization/stats
# pass over data the Python notebook already generated.
#
# Run from the same working directory as the notebook's CSV output, e.g.:
#   Rscript sensitivity_analysis.R

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

poi_names <- c("Rhododendron", "Parkdale", "Govt. Camp")

param_labels <- c(
  plume_height   = "Plume height (m asl)",
  eruption_mass  = "Eruption mass (kg)",
  diffusion_coef = "Diffusion coefficient (m^2/s)"
)

# Hazard thresholds (kg/m^2) -- Wilson et al. (2014), Jenkins et al. (2015). Drawn as reference
# lines on the OAT and tornado plots, matching tephra2_sensitivity.ipynb.
hazard_thresholds <- c(1, 10, 100, 1000)

# ---------------------------------------------------------------------------
# Load the three OAT sweeps
# ---------------------------------------------------------------------------
oat_files <- c(
  plume_height   = "oat_plume_height.csv",
  eruption_mass  = "oat_eruption_mass.csv",
  diffusion_coef = "oat_diffusion_coef.csv"
)

oat_data <- lapply(oat_files, read_csv, show_col_types = FALSE)

# ---------------------------------------------------------------------------
# Per-parameter plot: one line per POI, log-scale y (matches the notebook's
# plot_oat() -- x-axis log-scaled too if the parameter's own values span
# more than one order of magnitude, since that's how they'd have been swept)
# ---------------------------------------------------------------------------
plot_oat <- function(param_name) {
  df <- oat_data[[param_name]]
  df_long <- df %>%
    pivot_longer(cols = all_of(poi_names), names_to = "location", values_to = "mass_kg_m2")

  x_range_ratio <- max(df[[param_name]]) / min(df[[param_name]])
  use_log_x <- x_range_ratio > 10

  p <- ggplot(df_long, aes(x = .data[[param_name]], y = mass_kg_m2, color = location)) +
    geom_hline(yintercept = hazard_thresholds, color = "gray50", linetype = "dotted", linewidth = 0.4) +
    annotate("text", x = max(df[[param_name]]), y = hazard_thresholds,
             label = paste0(hazard_thresholds, " kg/m^2"), hjust = 1, vjust = -0.3,
             size = 2.5, color = "gray40") +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.2) +
    scale_y_log10() +
    labs(
      x = param_labels[[param_name]],
      y = "Ash thickness (kg/m^2)",
      color = "Location",
      title = paste0("OAT sensitivity: ", param_labels[[param_name]]),
      subtitle = "Other two parameters held at baseline"
    ) +
    theme_minimal()

  if (use_log_x) p <- p + scale_x_log10()

  ggsave(paste0("R_oat_", param_name, ".png"), p, width = 8, height = 5, dpi = 150)
  p
}

oat_plots <- lapply(names(oat_data), plot_oat)

# ---------------------------------------------------------------------------
# Summary: range (kg/m^2) and max/min ratio, per parameter PER POI, plus a
# sensitivity ranking at each location (not just one blended average) --
# matches the spec's requested summary shape.
# ---------------------------------------------------------------------------
summarize_oat <- function(param_name) {
  oat_data[[param_name]] %>%
    pivot_longer(cols = all_of(poi_names), names_to = "location", values_to = "mass_kg_m2") %>%
    group_by(location) %>%
    summarise(
      parameter = param_labels[[param_name]],
      min_mass_kg_m2 = min(mass_kg_m2, na.rm = TRUE),
      max_mass_kg_m2 = max(mass_kg_m2, na.rm = TRUE),
      range_kg_m2 = max_mass_kg_m2 - min_mass_kg_m2,
      max_min_ratio = ifelse(min_mass_kg_m2 > 0, max_mass_kg_m2 / min_mass_kg_m2, NA_real_),
      .groups = "drop"
    )
}

sensitivity_summary <- bind_rows(lapply(names(oat_data), summarize_oat)) %>%
  group_by(location) %>%
  mutate(rank_at_location = rank(-max_min_ratio, ties.method = "min")) %>%
  ungroup() %>%
  arrange(location, rank_at_location) %>%
  select(location, parameter, min_mass_kg_m2, max_mass_kg_m2, range_kg_m2,
         max_min_ratio, rank_at_location)

write_csv(sensitivity_summary, "R_sensitivity_summary.csv")
print(sensitivity_summary)

# ---------------------------------------------------------------------------
# Tornado diagram: one bar per parameter, spanning its min-max output range
# across all 3 POIs combined, sorted by range (widest = most sensitive).
# ---------------------------------------------------------------------------
tornado_data <- sensitivity_summary %>%
  group_by(parameter) %>%
  summarise(
    min_mass_kg_m2 = min(min_mass_kg_m2),
    max_mass_kg_m2 = max(max_mass_kg_m2),
    .groups = "drop"
  ) %>%
  arrange(max_mass_kg_m2 / pmax(min_mass_kg_m2, 1e-12))

tornado_data$parameter <- factor(tornado_data$parameter, levels = tornado_data$parameter)

p_tornado <- ggplot(tornado_data, aes(y = parameter)) +
  geom_vline(xintercept = hazard_thresholds, color = "gray50", linetype = "dotted", linewidth = 0.4) +
  geom_segment(aes(x = min_mass_kg_m2, xend = max_mass_kg_m2, yend = parameter),
               linewidth = 10, color = "#d95f02", alpha = 0.75, lineend = "butt") +
  scale_x_log10() +
  labs(
    x = "Ash thickness range across swept parameter (kg/m^2), all POIs combined",
    y = NULL,
    title = "Tornado diagram: which parameter moves ash thickness the most?"
  ) +
  theme_minimal()

ggsave("R_sensitivity_tornado.png", p_tornado, width = 8, height = 4, dpi = 150)

cat("\nWrote: R_oat_plume_height.png, R_oat_eruption_mass.png, R_oat_diffusion_coef.png,\n")
cat("       R_sensitivity_tornado.png, R_sensitivity_summary.csv\n")
