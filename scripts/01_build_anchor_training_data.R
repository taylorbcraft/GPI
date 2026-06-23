# Build the polygon-level anchor training table.
# Join field measurements, plant richness, and raster summaries for sampled polygons.

source("config.R")

library(tidyverse)
library(sf)
library(terra)
library(exactextractr)
library(janitor)

normalize_gpi_label <- function(x) {
  x %>%
    stringr::str_trim() %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[- ]+", "_") %>%
    stringr::str_replace("^mid_low$", "extensive") %>%
    stringr::str_replace("^mid_high$", "mid")
}

ensure_dirs(paths$training_dir)

# clean observer labels
env <- read_csv(paths$raw_env, show_col_types = FALSE) %>%
  clean_names() %>%
  rename(observer_estimated_GPI = in_lui) %>%
  mutate(
    polygon_id = as.character(polygon_id),
    observer_estimated_GPI = normalize_gpi_label(observer_estimated_GPI),
    observer_estimated_GPI = factor(observer_estimated_GPI, levels = gpi_class_levels)
  )

plant_div <- read_csv(paths$raw_plant, show_col_types = FALSE) %>%
  clean_names()

# sampled zone geometry
zones <- st_read(paths$sampled_zone_geometry, quiet = TRUE) %>%
  clean_names() %>%
  mutate(
    polygon_id = as.character(polygon_id),
    meadow_id = str_remove(polygon_id, "_.*$")
  )

# calibration raster layers
predictor_rasters <- predictor_bands %>%
  set_names() %>%
  map(~ rast(path_calibration_raster(.x)))

if ("s2rep" %in% names(predictor_rasters)) {
  predictor_rasters$s2rep <- clamp(predictor_rasters$s2rep, lower = 600, upper = 850, values = TRUE)
}

zones <- st_transform(zones, crs(predictor_rasters[[1]]))

# polygon field summaries
env_summary <- env %>%
  group_by(polygon_id) %>%
  summarise(
    observer_estimated_GPI = dplyr::first(na.omit(as.character(observer_estimated_GPI)), default = NA_character_),
    soil_moisture = mean(sm_mean, na.rm = TRUE),
    soil_resistance = mean(resistance_mean, na.rm = TRUE),
    vegetation_height = mean(vh_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(observer_estimated_GPI = factor(observer_estimated_GPI, levels = gpi_class_levels))

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

write_csv(anchor_zone_training_data, path_anchor_training())
