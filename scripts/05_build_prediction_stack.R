# Build the prediction-year raster stack used by the pixel-level classifier.

source("config.R")

library(tidyverse)
library(terra)

ensure_dirs(paths$prediction_dir)

# Load only model bands
predictor_stack <- load_predictor_stack(
  band_names = model_predictor_bands,
  image_date = prediction_image_date
)

writeRaster(
  predictor_stack,
  filename = path_prediction_stack(),
  overwrite = TRUE
)
