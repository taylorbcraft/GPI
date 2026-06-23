# Run the full calibration and prediction pipeline.
# This rebuilds the training products, fits the model, and writes the final map outputs.

scripts_to_run <- c(
  "scripts/01_build_anchor_training_data.R",
  "scripts/02_validate_environmental_relationships.R",
  "scripts/03_define_candidate_gpi_classes.R",
  "scripts/04_train_gpi_classifier.R",
  "scripts/05_predict_pixel_map.R"
)

for (script_path in scripts_to_run) {
  message("running ", script_path)
  source(script_path)
}
