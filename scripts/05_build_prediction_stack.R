# Build the prediction-year raster stack used by the pixel-level classifier.

source("config.R")

library(tidyverse)
library(terra)

ensure_dirs(paths$prediction_dir)

# Load only model bands
predictor_stack <- model_predictor_bands %>%
  set_names() %>%
  map(~ rast(path_prediction_raster(.x)))

if ("s2rep" %in% names(predictor_stack)) {
  predictor_stack$s2rep <- clamp(predictor_stack$s2rep, lower = 600, upper = 850, values = TRUE)
}

predictor_stack <- rast(predictor_stack)
names(predictor_stack) <- model_predictor_bands

writeRaster(
  predictor_stack,
  filename = path_prediction_stack(),
  overwrite = TRUE
)
