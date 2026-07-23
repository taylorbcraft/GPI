# Run the full calibration and prediction pipeline.
# The scripts are sourced in dependency order to rebuild training data, model
# validation, the field habitat map, and the analysis figures and tables.

# run workflow
for (script_path in c(
  "scripts/01_build_anchor_training_data.R",
  "scripts/02_validate_environmental_relationships.R",
  "scripts/03_define_candidate_habitat_classes.R",
  "scripts/04_train_habitat_classifier.R",
  "scripts/05_predict_habitat_map.R",
  "scripts/06_compare_s1_s2_classes.R",
  "scripts/07_analyze_biotic_relationships.R"
)) {
  message("running ", script_path)
  source(script_path)
}
