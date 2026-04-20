# Build the field-level predictor table used for final mapping.
# Extracts the model raster means for every polygon in the full field layer.
#
# Inputs: field_geometry.gpkg and processed rasters.
# Output: data/processed/predictions/field_predictor_data_<prediction_year>.csv

source("config.R")

library(tidyverse)
library(sf)
library(terra)
library(exactextractr)
library(janitor)

# Load only the prediction-year bands needed by the saved model.
load_predictor_rasters <- function(band_names) {
  rasters <- band_names %>%
    set_names() %>%
    map(~ rast(path_prediction_raster(.x)))

  if ("s2rep" %in% names(rasters)) {
    # Clip S2REP formula blow-ups before polygon summaries.
    # The 600-850 nm bounds are broad guards around plausible red-edge positions.
    rasters$s2rep <- clamp(rasters$s2rep, lower = 600, upper = 850, values = TRUE)
  }

  rasters
}

ensure_dirs(paths$prediction_dir)

require_file(paths$field_geometry)
purrr::walk(model_predictor_bands, ~ require_file(path_prediction_raster(.x), paste("prediction raster for", .x)))

# The field layer is the mapping unit; every retained field receives predictors.
fields <- st_read(paths$field_geometry, quiet = TRUE) %>%
  clean_names()

if (!field_id_col %in% names(fields)) {
  stop(paste0("Field identifier column not found in field geometry: ", field_id_col), call. = FALSE)
}

names(fields)[names(fields) == field_id_col] <- "field_id"

fields <- fields %>%
  mutate(field_id = as.character(field_id))

predictor_rasters <- load_predictor_rasters(model_predictor_bands)

# Match the field geometry to the raster CRS before extraction.
fields <- st_transform(fields, crs(predictor_rasters[[1]]))

extract_mean <- function(r) {
  exact_extract(r, fields, "mean")
}

# Summarize prediction-year rasters using the same model predictors as training.
predictor_summary <- imap_dfc(
  predictor_rasters,
  ~ tibble(!!.y := extract_mean(.x))
)

field_predictor_data <- fields %>%
  st_drop_geometry() %>%
  select(field_id) %>%
  bind_cols(predictor_summary)

write_csv(field_predictor_data, path_field_predictor_data())
