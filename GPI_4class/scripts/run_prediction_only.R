# Rebuild the prediction stack and apply the saved model to the current rasters.

scripts_to_run <- c(
  "scripts/05_build_prediction_stack.R",
  "scripts/06_predict_pixel_map.R"
)

for (script_path in scripts_to_run) {
  message("running ", script_path)
  source(script_path)
}
