# Run the prediction-only pipeline.
# This rebuilds the prediction stack and applies the saved model to current rasters.

scripts_to_run <- c(
  "scripts/05_predict_pixel_map.R"
)

for (script_path in scripts_to_run) {
  message("running ", script_path)
  source(script_path)
}
