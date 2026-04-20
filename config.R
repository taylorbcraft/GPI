# Shared settings for the GPI pipeline.
# The 2025 field data calibrate the model. Change prediction_year and
# prediction_image_date to map another year with the saved calibration model.

calibration_year <- "2025"
calibration_image_date <- "2025-04-11"

# Change these two values to make an annual map with the saved calibration model.
prediction_year <- "2025"
prediction_image_date <- "2025-04-11"

field_id_col <- "meadow_id"
zone_id_col <- "polygon_id"
meadow_id_col <- "meadow_id"
observer_gpi_col <- "observer_estimated_GPI"

predictor_bands <- c(
  "s2rep",
  "ndvi",
  "ndwi",
  "savi",
  "evi",
  "msi",
  "ndmi",
  "mndwi"
)

# Only these bands enter the final random forest model.
model_predictor_bands <- c(
  "s2rep",
  "msi",
  "ndmi",
  "mndwi"
)

# KNN turns field measurements into the supervised training target.
knn_k <- 7
knn_vote_method <- "weighted_distance"

# Final mapped classes after collapsing the original four observer labels.
gpi_class_levels <- c(
  "extensive",
  "mid",
  "intensive"
)

paths <- list(
  raw_env = file.path("data", "raw", paste0("environmental_field_data_", calibration_year, ".csv")),
  raw_plant = file.path("data", "raw", paste0("plant_diversity_plots_", calibration_year, ".csv")),
  sampled_zone_geometry = file.path("data", "spatial", "sampled_zone_geometry.gpkg"),
  field_geometry = file.path("data", "spatial", "field_geometry.gpkg"),
  raster_dir = file.path("data", "processed", "rasters"),
  training_dir = file.path("data", "processed", "training"),
  validation_dir = file.path("data", "processed", "validation"),
  model_dir = file.path("data", "processed", "models"),
  prediction_dir = file.path("data", "processed", "predictions"),
  spatial_dir = file.path("data", "processed", "spatial"),
  figure_dir = "figures"
)

path_raster <- function(band_name, image_date) {
  file.path(paths$raster_dir, paste0(band_name, "_", image_date, "_mosaic.tif"))
}

path_calibration_raster <- function(band_name) {
  path_raster(band_name, calibration_image_date)
}

path_prediction_raster <- function(band_name) {
  path_raster(band_name, prediction_image_date)
}

path_anchor_training <- function() {
  file.path(paths$training_dir, paste0("anchor_zone_training_data_", calibration_year, ".csv"))
}

path_candidate_training <- function() {
  file.path(paths$training_dir, paste0("candidate_gpi_training_data_", calibration_year, ".csv"))
}

path_environmental_validation_summary <- function() {
  file.path(paths$validation_dir, paste0("environmental_validation_summary_", calibration_year, ".csv"))
}

path_environmental_validation_plot <- function() {
  file.path(paths$figure_dir, paste0("environmental_validation_plots_", calibration_year, ".png"))
}

path_candidate_class_plot <- function() {
  file.path(paths$figure_dir, paste0("candidate_gpi_class_boxplots_", calibration_year, ".png"))
}

path_estimated_thresholds <- function() {
  file.path(paths$validation_dir, paste0("gpi_estimated_rule_thresholds_", calibration_year, ".csv"))
}

path_model_comparison <- function() {
  file.path(paths$validation_dir, paste0("gpi_model_comparison_", calibration_year, ".csv"))
}

path_model_tuning <- function() {
  file.path(paths$validation_dir, paste0("gpi_model_tuning_", calibration_year, ".csv"))
}

path_best_model_confusion <- function() {
  file.path(paths$validation_dir, paste0("gpi_best_model_confusion_matrix_", calibration_year, ".csv"))
}

path_best_model_class_accuracy <- function() {
  file.path(paths$validation_dir, paste0("gpi_best_model_class_accuracy_", calibration_year, ".csv"))
}

path_best_model_variable_importance <- function() {
  file.path(paths$validation_dir, paste0("gpi_best_model_variable_importance_", calibration_year, ".csv"))
}

path_best_model_metadata <- function() {
  file.path(paths$validation_dir, paste0("gpi_best_model_metadata_", calibration_year, ".csv"))
}

path_best_model <- function() {
  file.path(paths$model_dir, paste0("gpi_best_model_", calibration_year, ".rds"))
}

path_field_predictor_data <- function() {
  file.path(paths$prediction_dir, paste0("field_predictor_data_", prediction_year, ".csv"))
}

path_field_predictions <- function() {
  file.path(paths$prediction_dir, paste0("field_gpi_predictions_", prediction_year, ".csv"))
}

path_field_map <- function() {
  file.path(paths$spatial_dir, paste0("field_gpi_map_", prediction_year, ".gpkg"))
}

path_field_map_preview <- function() {
  file.path(paths$figure_dir, paste0("field_gpi_map_", prediction_year, ".png"))
}

ensure_dirs <- function(dir_paths) {
  invisible(lapply(dir_paths, dir.create, recursive = TRUE, showWarnings = FALSE))
}

require_file <- function(path, label = NULL) {
  if (!file.exists(path)) {
    if (is.null(label)) {
      label <- path
    }
    stop(paste0("Missing required file: ", label, " (", path, ")"), call. = FALSE)
  }
}
