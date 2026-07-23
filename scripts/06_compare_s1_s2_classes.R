# Compare Sentinel-1 temporal variability with field-level Sentinel-2 habitat classes.
# The script extracts mean Sentinel-1 values to each field polygon, summarizes
# those values by dominant Sentinel-2 habitat class, runs rank-based tests, and saves
# analysis outputs plus a comparison figure.

source("config.R")

library(tidyverse)
library(sf)
library(terra)
library(exactextractr)
library(ggplot2)

# settings
s1_year <- Sys.getenv("S1_YEAR", unset = prediction_year)
s1_polarization <- Sys.getenv("S1_POLARIZATION", unset = "VV")
s1_orbit_pass <- Sys.getenv("S1_ORBIT_PASS", unset = "ASCENDING")

s1_raster_path <- file.path(
  paths$raster_dir,
  paste0(
    "S1_",
    s1_polarization,
    "_LogRatio_StdDev_",
    str_to_title(str_to_lower(s1_orbit_pass)),
    "_",
    s1_year,
    ".tif"
  )
)

ensure_dirs(c(paths$validation_dir, paths$figure_dir))

# load data
field_habitat_map <- st_read(path_field_map(prediction_year), quiet = TRUE) %>%
  mutate(
    habitat_class = factor(habitat_class, levels = habitat_class_levels),
    class_rank = match(habitat_class, habitat_class_levels)
  )

s1_raster <- rast(s1_raster_path)

if (!identical(st_crs(field_habitat_map)$wkt, crs(s1_raster))) {
  field_habitat_map <- st_transform(field_habitat_map, crs(s1_raster))
}

# extract field means
field_values <- exact_extract(
  s1_raster,
  field_habitat_map,
  "mean",
  progress = FALSE
) %>%
  tibble(s1_temporal_sd = .) %>%
  bind_cols(field_habitat_map %>% st_drop_geometry() %>% select(field_id, habitat_class, class_rank)) %>%
  relocate(field_id, habitat_class, class_rank, s1_temporal_sd) %>%
  filter(!is.na(habitat_class), !is.na(s1_temporal_sd), !is.na(class_rank))

# summarize classes
class_summary <- field_values %>%
  group_by(habitat_class, class_rank) %>%
  summarise(
    n_fields = n(),
    mean_s1_temporal_sd = mean(s1_temporal_sd),
    sd_s1_temporal_sd = sd(s1_temporal_sd),
    median_s1_temporal_sd = median(s1_temporal_sd),
    iqr_s1_temporal_sd = IQR(s1_temporal_sd),
    min_s1_temporal_sd = min(s1_temporal_sd),
    max_s1_temporal_sd = max(s1_temporal_sd),
    .groups = "drop"
  ) %>%
  arrange(class_rank) %>%
  select(-class_rank)

# tests
kruskal_result <- kruskal.test(s1_temporal_sd ~ habitat_class, data = field_values)
pairwise_result <- pairwise.wilcox.test(
  x = field_values$s1_temporal_sd,
  g = field_values$habitat_class,
  p.adjust.method = "holm"
)
spearman_result <- cor.test(
  ~ s1_temporal_sd + class_rank,
  data = field_values,
  method = "spearman"
)

# tabular test output
pairwise_tbl <- as_tibble(as.table(pairwise_result$p.value), .name_repair = "minimal") %>%
  set_names(c("class_1", "class_2", "p_value")) %>%
  filter(!is.na(p_value)) %>%
  mutate(
    test = "pairwise_wilcoxon",
    comparison = paste(class_1, "vs", class_2),
    statistic = NA_real_,
    df = NA_real_,
    estimate = NA_real_,
    estimate_label = NA_character_,
    p_adjustment = "holm",
    .before = 1
  ) %>%
  select(test, comparison, statistic, df, estimate, estimate_label, p_value, p_adjustment)

test_summary <- bind_rows(
  tibble(
    test = "kruskal_wallis",
    comparison = "all classes",
    statistic = unname(kruskal_result$statistic),
    df = unname(kruskal_result$parameter),
    estimate = NA_real_,
    estimate_label = NA_character_,
    p_value = kruskal_result$p.value,
    p_adjustment = NA_character_
  ),
  pairwise_tbl,
  tibble(
    test = "spearman_rank_correlation",
    comparison = "class rank vs Sentinel-1 temporal SD",
    statistic = unname(spearman_result$statistic),
    df = NA_real_,
    estimate = unname(spearman_result$estimate),
    estimate_label = "rho",
    p_value = spearman_result$p.value,
    p_adjustment = NA_character_
  )
)

# save analysis outputs
write_csv(
  field_values,
  file.path(paths$validation_dir, paste0("s1_s2_field_values_", s1_polarization, "_", s1_year, ".csv"))
)
write_csv(
  class_summary,
  file.path(paths$validation_dir, paste0("s1_s2_class_comparison_summary_", s1_polarization, "_", s1_year, ".csv"))
)
write_csv(
  test_summary,
  file.path(paths$validation_dir, paste0("s1_s2_class_comparison_tests_", s1_polarization, "_", s1_year, ".csv"))
)

stats_lines <- c(
  paste("Sentinel-1 raster:", s1_raster_path),
  paste("Field map:", path_field_map(prediction_year)),
  paste("Fields compared:", nrow(field_values)),
  "",
  "Kruskal-Wallis test",
  paste("chi-squared =", unname(kruskal_result$statistic)),
  paste("df =", unname(kruskal_result$parameter)),
  paste("p-value =", kruskal_result$p.value),
  "",
  "Pairwise Wilcoxon tests (Holm-adjusted p-values)",
  capture.output(print(pairwise_result$p.value)),
  "",
  "Spearman rank correlation",
  paste("rho =", unname(spearman_result$estimate)),
  paste("S =", unname(spearman_result$statistic)),
  paste("p-value =", spearman_result$p.value)
)

write_lines(
  stats_lines,
  file.path(paths$validation_dir, paste0("s1_s2_class_comparison_stats_", s1_polarization, "_", s1_year, ".txt"))
)

# save comparison plot
comparison_plot <- ggplot(
  field_values,
  aes(x = habitat_class, y = s1_temporal_sd, fill = habitat_class)) +
  coord_cartesian(ylim = c(0.1, 0.25)) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.85,
    staplewidth = 0.5,
    linewidth = 0.8
  ) +
  scale_fill_manual(
    values = habitat_class_palette,
    breaks = habitat_class_levels,
    drop = FALSE
  ) +
  labs(
    x = "Habitat quality",
    y = "Sentinel-1 LUI",
    fill = "Habitat quality"
  ) +
  theme_minimal(base_family = "Avenir Next", base_size = 12) +
  theme(
    text = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.4),
    axis.text = element_text(color = "grey35"),
    axis.title = element_text(color = "grey20")
  )

ggsave(
  filename = file.path(paths$figure_dir, paste0("s1_s2_class_comparison_", s1_polarization, "_", s1_year, ".png")),
  plot = comparison_plot,
  device = ragg::agg_png,
  width = 7,
  height = 5,
  dpi = 300
)
