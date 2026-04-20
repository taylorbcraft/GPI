# Apply the saved calibration model to the configured prediction year.
# Requires prediction rasters and field geometry, but not new field observations.

scripts_to_run <- c(
  "scripts/05_build_field_prediction_data.R",
  "scripts/06_predict_field_gpi_classes.R"
)

for (script_path in scripts_to_run) {
  message("running ", script_path)
  source(script_path)
}
