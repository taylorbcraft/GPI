# Compare Sentinel-1 temporal variability with field-level Sentinel-2 classes.
# The script extracts mean Sentinel-1 values to each field polygon, summarizes
# those values by dominant Sentinel-2 class, runs rank-based tests, and saves
# table outputs plus a comparison figure.

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
field_gpi_map <- st_read(path_field_map(), quiet = TRUE) %>%
  mutate(
    dominant_class = factor(dominant_class, levels = gpi_class_levels),
    class_rank = match(dominant_class, gpi_class_levels)
  )

s1_raster <- rast(s1_raster_path)

if (!identical(st_crs(field_gpi_map)$wkt, crs(s1_raster))) {
  field_gpi_map <- st_transform(field_gpi_map, crs(s1_raster))
}

# extract field means
field_values <- exact_extract(
  s1_raster,
  field_gpi_map,
  "mean",
  progress = FALSE
) %>%
  tibble(s1_temporal_sd = .) %>%
  bind_cols(field_gpi_map %>% st_drop_geometry() %>% select(field_id, dominant_class, class_rank)) %>%
  relocate(field_id, dominant_class, class_rank, s1_temporal_sd) %>%
  filter(!is.na(dominant_class), !is.na(s1_temporal_sd), !is.na(class_rank))

# summarize classes
class_summary <- field_values %>%
  group_by(dominant_class, class_rank) %>%
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
kruskal_result <- kruskal.test(s1_temporal_sd ~ dominant_class, data = field_values)
pairwise_result <- pairwise.wilcox.test(
  x = field_values$s1_temporal_sd,
  g = field_values$dominant_class,
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

# write tables and stats
write_csv(field_values, path_s1_s2_field_values(s1_polarization, s1_year))
write_csv(class_summary, path_s1_s2_class_comparison_summary(s1_polarization, s1_year))
write_csv(test_summary, path_s1_s2_class_comparison_tests(s1_polarization, s1_year))

stats_lines <- c(
  paste("Sentinel-1 raster:", s1_raster_path),
  paste("Field map:", path_field_map()),
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
  path_s1_s2_class_comparison_stats(s1_polarization, s1_year)
)

# save comparison plot
comparison_plot <- ggplot(
  field_values,
  aes(x = dominant_class, y = s1_temporal_sd, fill = dominant_class)
) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.18, height = 0, size = 0.6, alpha = 0.2, color = "black") +
  scale_fill_manual(
    values = c(extensive = "#2E7D32", mid = "#FDAE61", intensive = "#D73027"),
    breaks = gpi_class_levels,
    drop = FALSE
  ) +
  labs(
    x = "Sentinel-2 habitat quality class",
    y = "Sentinel-1 LUI",
    fill = "Class"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

comparison_plot

ggsave(
  filename = path_s1_s2_class_comparison_plot(s1_polarization, s1_year),
  plot = comparison_plot,
  width = 7,
  height = 5,
  dpi = 300
)
