# Rebuild the 2025 calibration model and map the configured prediction year.
# Use this after calibration data, rasters, or model settings change.

scripts_to_run <- c(
  "scripts/01_build_anchor_training_data.R",
  "scripts/02_validate_environmental_relationships.R",
  "scripts/03_define_candidate_gpi_classes.R",
  "scripts/04_train_gpi_classifier.R",
  "scripts/05_build_field_prediction_data.R",
  "scripts/06_predict_field_gpi_classes.R"
)

for (script_path in scripts_to_run) {
  message("running ", script_path)
  source(script_path)
}
