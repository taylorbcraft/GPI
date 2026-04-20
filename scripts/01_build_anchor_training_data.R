# Build the anchor training table.
# Joins sampled-zone geometry, field observations, plant richness, and raster
# means at polygon_id. Script 02 validates these relationships, and script 03
# turns the observer labels into the KNN-derived GPI target.
#
# Inputs: raw field CSVs, sampled_zone_geometry.gpkg, and processed rasters.
# Output: data/processed/training/anchor_zone_training_data_<calibration_year>.csv

source("config.R")

library(tidyverse)
library(sf)
library(terra)
library(exactextractr)
library(janitor)

# Standardize observer labels and apply the final three-class scheme at intake.
# Downstream scripts should only see extensive, mid, and intensive.
normalize_gpi_label <- function(x) {
  x %>%
    stringr::str_trim() %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[- ]+", "_") %>%
    stringr::str_replace("^mid_low$", "extensive") %>%
    stringr::str_replace("^mid_high$", "mid")
}

first_non_missing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_character_)
  }
  as.character(x[[1]])
}

# Load the calibration-year raster stack used to summarize sampled zones.
load_predictor_rasters <- function() {
  rasters <- predictor_bands %>%
    set_names() %>%
    map(~ rast(path_calibration_raster(.x)))

  if ("s2rep" %in% names(rasters)) {
    # Clip S2REP formula blow-ups before polygon summaries.
    # The 600-850 nm bounds are broad guards around plausible red-edge positions.
    rasters$s2rep <- clamp(rasters$s2rep, lower = 600, upper = 850, values = TRUE)
  }

  rasters
}

ensure_dirs(paths$training_dir)

require_file(paths$raw_env)
require_file(paths$raw_plant)
require_file(paths$sampled_zone_geometry)
purrr::walk(predictor_bands, ~ require_file(path_calibration_raster(.x), paste("calibration raster for", .x)))

# The source CSV still uses in_lui; rename it and collapse labels once here.
env <- read_csv(paths$raw_env, show_col_types = FALSE) %>%
  clean_names() %>%
  rename(observer_estimated_GPI = in_lui) %>%
  mutate(
    polygon_id = as.character(polygon_id),
    observer_estimated_GPI = normalize_gpi_label(observer_estimated_GPI)
  )

unexpected_observer_labels <- env %>%
  filter(!is.na(observer_estimated_GPI), !observer_estimated_GPI %in% gpi_class_levels) %>%
  distinct(observer_estimated_GPI)

if (nrow(unexpected_observer_labels) > 0) {
  stop(
    paste0(
      "Unexpected observer_estimated_GPI labels after normalization: ",
      paste(unexpected_observer_labels$observer_estimated_GPI, collapse = ", ")
    ),
    call. = FALSE
  )
}

env <- env %>%
  mutate(
    observer_estimated_GPI = factor(observer_estimated_GPI, levels = gpi_class_levels)
  )

plant_div <- read_csv(paths$raw_plant, show_col_types = FALSE) %>%
  clean_names()

# Sampled-zone geometry defines the polygon grain for field and raster summaries.
zones <- st_read(paths$sampled_zone_geometry, quiet = TRUE) %>%
  clean_names() %>%
  mutate(
    polygon_id = as.character(polygon_id),
    meadow_id = str_remove(polygon_id, "_.*$")
  )

# Each sampled zone should carry at most one observer class.
observer_conflicts <- env %>%
  filter(!is.na(observer_estimated_GPI)) %>%
  group_by(polygon_id) %>%
  summarise(n_labels = n_distinct(observer_estimated_GPI), .groups = "drop") %>%
  filter(n_labels > 1)

if (nrow(observer_conflicts) > 0) {
  stop(
    paste0(
      "Conflicting observer_estimated_GPI labels found within polygon_id values: ",
      paste(observer_conflicts$polygon_id, collapse = ", ")
    ),
    call. = FALSE
  )
}

predictor_rasters <- load_predictor_rasters()

# Align sampled zones to the raster CRS before extraction.
zones <- st_transform(zones, crs(predictor_rasters[[1]]))

# Collapse repeated field measurements to one record per sampled polygon.
env_summary <- env %>%
  group_by(polygon_id) %>%
  summarise(
    observer_estimated_GPI = first_non_missing(observer_estimated_GPI),
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

extract_mean <- function(r) {
  exact_extract(r, zones, "mean")
}

# Extract raster means at the same polygon grain as the field observations.
predictor_summary <- imap_dfc(
  predictor_rasters,
  ~ tibble(!!.y := extract_mean(.x))
)

rs_summary <- zones %>%
  st_drop_geometry() %>%
  select(polygon_id, meadow_id) %>%
  bind_cols(predictor_summary)

# Join ecological and raster summaries into the anchor table used downstream.
anchor_zone_training_data <- rs_summary %>%
  left_join(env_summary, by = "polygon_id") %>%
  left_join(plant_summary, by = "polygon_id")

write_csv(anchor_zone_training_data, path_anchor_training())
