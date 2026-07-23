# Build the polygon-level anchor training table.
# The script standardizes observer labels, aggregates soil and vegetation surveys,
# extracts April Sentinel-2 predictor means for each sampled polygon, and joins
# those inputs into the calibration dataset.

source("config.R")

library(tidyverse)
library(sf)
library(terra)
library(exactextractr)
library(janitor)

normalize_habitat_class_label <- function(x) {
  x %>%
    stringr::str_trim() %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[- ]+", "_") %>%
    stringr::str_replace("^mid_low$", "high") %>%
    stringr::str_replace("^mid_high$", "moderate") %>%
    stringr::str_replace("^extensive$", "high") %>%
    stringr::str_replace("^mid$", "moderate") %>%
    stringr::str_replace("^intensive$", "low")
}

ensure_dirs(paths$training_dir)

# clean observer labels
env <- read_csv(
  file.path("data", "raw", "environmental_field_data_2025.csv"),
  show_col_types = FALSE
) %>%
  clean_names() %>%
  rename(observer_estimated_class = in_lui) %>%
  mutate(
    polygon_id = as.character(polygon_id),
    observer_estimated_class = normalize_habitat_class_label(observer_estimated_class),
    observer_estimated_class = factor(observer_estimated_class, levels = habitat_class_levels)
  )

plant_div <- read_csv(
  file.path("data", "raw", "plant_diversity_plots_2025.csv"),
  show_col_types = FALSE
) %>%
  clean_names()

# sampled zone geometry
zones <- st_read(paths$sampled_zone_geometry, quiet = TRUE) %>%
  clean_names() %>%
  mutate(
    polygon_id = as.character(polygon_id),
    meadow_id = str_remove(polygon_id, "_.*$")
  )

# calibration raster layers
predictor_rasters <- candidate_model_predictor_bands %>%
  set_names() %>%
  map(~ rast(path_raster(.x, "2025-04-median")))

if ("s2rep" %in% names(predictor_rasters)) {
  predictor_rasters$s2rep <- clamp(predictor_rasters$s2rep, lower = 600, upper = 850, values = TRUE)
}

zones <- st_transform(zones, crs(predictor_rasters[[1]]))

# polygon field summaries
env_summary <- env %>%
  group_by(polygon_id) %>%
  summarise(
    observer_estimated_class = dplyr::first(na.omit(as.character(observer_estimated_class)), default = NA_character_),
    soil_moisture = mean(sm_mean, na.rm = TRUE),
    soil_resistance = mean(resistance_mean, na.rm = TRUE),
    vegetation_height = mean(vh_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(observer_estimated_class = factor(observer_estimated_class, levels = habitat_class_levels))

plant_summary <- plant_div %>%
  group_by(polygon_id) %>%
  summarise(
    plant_richness_mean = mean(spec_count_plot, na.rm = TRUE),
    .groups = "drop"
  )

# polygon raster means
predictor_summary <- imap_dfc(
  predictor_rasters,
  ~ tibble(!!.y := exact_extract(.x, zones, "mean"))
)

anchor_zone_training_data <- zones %>%
  st_drop_geometry() %>%
  select(polygon_id, meadow_id) %>%
  bind_cols(predictor_summary) %>%
  left_join(env_summary, by = "polygon_id") %>%
  left_join(plant_summary, by = "polygon_id")

# save training data
write_csv(anchor_zone_training_data, path_anchor_training())
