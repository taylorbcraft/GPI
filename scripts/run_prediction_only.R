# Run the prediction-only pipeline.
# The script applies the saved model to the configured prediction-year rasters
# and rebuilds the field-level habitat map and map figure.

# run prediction
source("scripts/05_predict_habitat_map.R")
