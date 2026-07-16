# Shared settings for the GPI pipeline.
# Update the calibration and prediction dates here before running the scripts.

# calibration year and imagery
# default workflow uses the April monthly median mosaic
calibration_year <- "2025"
calibration_image_date <- "2025-04-median"

# prediction year and imagery
prediction_year <- "2025"
prediction_image_date <- "2025-04-median"

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

# candidate predictor bands considered in rf selection
candidate_model_predictor_bands <- predictor_bands

# max predictors retained by forward selection
max_model_predictors <- 2

# predictor aggregation unit used in classifier training
model_training_unit <- "pixel"

# final ordered class labels
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

custom_predictor_rasters <- c()

path_raster <- function(band_name, image_date) {
  if (band_name %in% names(custom_predictor_rasters)) {
    return(custom_predictor_rasters[[band_name]])
  }

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

path_environmental_validation_coefficients <- function() {
  file.path(paths$validation_dir, paste0("environmental_validation_coefficients_", calibration_year, ".csv"))
}

path_environmental_validation_full_summary <- function() {
  file.path(paths$validation_dir, paste0("environmental_validation_full_summary_", calibration_year, ".txt"))
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

path_s1_s2_class_comparison_summary <- function(polarization = "VV", year = prediction_year) {
  file.path(
    paths$validation_dir,
    paste0("s1_s2_class_comparison_summary_", polarization, "_", year, ".csv")
  )
}

path_s1_s2_class_comparison_stats <- function(polarization = "VV", year = prediction_year) {
  file.path(
    paths$validation_dir,
    paste0("s1_s2_class_comparison_stats_", polarization, "_", year, ".txt")
  )
}

path_s1_s2_class_comparison_tests <- function(polarization = "VV", year = prediction_year) {
  file.path(
    paths$validation_dir,
    paste0("s1_s2_class_comparison_tests_", polarization, "_", year, ".csv")
  )
}

path_s1_s2_field_values <- function(polarization = "VV", year = prediction_year) {
  file.path(
    paths$validation_dir,
    paste0("s1_s2_field_values_", polarization, "_", year, ".csv")
  )
}

path_s1_s2_class_comparison_plot <- function(polarization = "VV", year = prediction_year) {
  file.path(
    paths$figure_dir,
    paste0("s1_s2_class_comparison_", polarization, "_", year, ".png")
  )
}

ensure_dirs <- function(dir_paths) {
  invisible(lapply(dir_paths, dir.create, recursive = TRUE, showWarnings = FALSE))
}

load_predictor_stack <- function(band_names, image_date, resample_method = "near") {
  raster_list <- vector("list", length(band_names))
  names(raster_list) <- band_names
  template <- NULL

  for (band_name in band_names) {
    raster_layer <- rast(path_raster(band_name, image_date))
    names(raster_layer) <- band_name

    if (is.null(template)) {
      template <- raster_layer
    } else if (!compareGeom(template, raster_layer, stopOnError = FALSE)) {
      if (!same.crs(template, raster_layer)) {
        raster_layer <- project(raster_layer, template, method = resample_method)
      } else {
        raster_layer <- resample(raster_layer, template, method = resample_method)
      }
    }

    if (band_name == "s2rep") {
      raster_layer <- clamp(raster_layer, lower = 600, upper = 850, values = TRUE)
    }

    raster_list[[band_name]] <- raster_layer
  }

  rast(raster_list)
}
