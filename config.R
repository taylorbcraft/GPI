# Define shared habitat-quality workflow settings.
# The file centralizes project paths, predictor names, class labels, plotting
# colors, and reusable path helpers sourced by the analysis scripts.

# project paths
paths <- list(
  sampled_zone_geometry = file.path("data", "spatial", "sampled_zone_geometry.gpkg"),
  field_geometry = file.path("data", "spatial", "field_geometry.gpkg"),
  raster_dir = file.path("data", "processed", "rasters"),
  training_dir = file.path("data", "processed", "training"),
  validation_dir = file.path("data", "processed", "validation"),
  model_dir = file.path("data", "processed", "models"),
  spatial_dir = file.path("data", "processed", "spatial"),
  figure_dir = "figures"
)

prediction_year <- "2025"

candidate_model_predictor_bands <- c(
  "s2rep",
  "ndvi",
  "ndwi",
  "savi",
  "evi",
  "msi",
  "ndmi",
  "mndwi"
)

habitat_class_levels <- c(
  "low",
  "moderate",
  "high"
)

habitat_class_palette <- c(
  low = "#CC79A7",
  moderate = "#F0E442",
  high = "#009E73"
)

# raster path
path_raster <- function(band_name, image_date) {
  file.path(paths$raster_dir, paste0(band_name, "_", image_date, "_mosaic.tif"))
}

# calibration paths
path_anchor_training <- function() {
  file.path(paths$training_dir, "anchor_zone_training_data_2025.csv")
}

path_candidate_training <- function() {
  file.path(paths$training_dir, "candidate_habitat_class_training_data_2025.csv")
}

path_best_model_metadata <- function() {
  file.path(paths$validation_dir, "habitat_best_model_metadata_2025.csv")
}

path_best_model <- function() {
  file.path(paths$model_dir, "habitat_best_model_2025.rds")
}

# field map path
path_field_map <- function(year) {
  file.path(paths$spatial_dir, paste0("field_habitat_class_map_", year, ".gpkg"))
}

# create directories
ensure_dirs <- function(dir_paths) {
  invisible(lapply(dir_paths, dir.create, recursive = TRUE, showWarnings = FALSE))
}
