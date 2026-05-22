# Shared settings for the 4-class GPI pipeline.
# This parallel project keeps the original observer classes separate instead of
# collapsing them into a 3-class system.

calibration_year <- "2025"
calibration_image_date <- "2025-04-11"

prediction_year <- "2025"
prediction_image_date <- "2025-04-11"

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

model_predictor_bands <- c(
  "s2rep",
  "msi",
  "ndmi",
  "mndwi"
)

knn_k <- 7

gpi_class_levels <- c(
  "extensive",
  "mid_low",
  "mid_high",
  "intensive"
)

gpi_class_display_labels <- c(
  extensive = "Extensive",
  mid_low = "Mid-low",
  mid_high = "Mid-high",
  intensive = "Intensive"
)

gpi_class_palette <- c(
  extensive = "#2E7D32",
  mid_low = "#A6D96A",
  mid_high = "#FDAE61",
  intensive = "#D73027"
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

path_prediction_stack <- function() {
  file.path(paths$prediction_dir, paste0("pixel_predictor_stack_", prediction_year, ".tif"))
}

path_pixel_predictions <- function() {
  file.path(paths$prediction_dir, paste0("field_gpi_pixel_summary_", prediction_year, ".csv"))
}

path_pixel_map <- function() {
  file.path(paths$spatial_dir, paste0("pixel_gpi_map_", prediction_year, ".tif"))
}

path_field_map <- function() {
  file.path(paths$spatial_dir, paste0("field_gpi_map_", prediction_year, ".gpkg"))
}

path_field_map_preview <- function() {
  file.path(paths$figure_dir, paste0("field_gpi_map_", prediction_year, ".png"))
}

path_pixel_map_preview <- function() {
  file.path(paths$figure_dir, paste0("pixel_gpi_map_", prediction_year, ".png"))
}

ensure_dirs <- function(dir_paths) {
  invisible(lapply(dir_paths, dir.create, recursive = TRUE, showWarnings = FALSE))
}
