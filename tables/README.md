# Tables

This folder contains copied tabular outputs that are likely to be cited,
summarized, or converted into tables in the manuscript.

Files

- `environmental_validation_summary_2025.csv`: summary fit metrics for the joint GAM relating S2REP to biophysical variables
- `environmental_validation_coefficients_2025.csv`: parametric and smooth-term output from the joint GAM
- `gpi_estimated_rule_thresholds_2025.csv`: class-definition metadata and thresholds for the rule-based classifier
- `gpi_model_comparison_2025.csv`: cross-validated performance of candidate random-forest target and predictor combinations
- `gpi_best_model_metadata_2025.csv`: metadata for the selected final random-forest model
- `gpi_best_model_confusion_matrix_2025.csv`: confusion matrix for the selected model
- `gpi_best_model_class_accuracy_2025.csv`: class-specific accuracy metrics for the selected model
- `gpi_best_model_variable_importance_2025.csv`: variable-importance scores for the selected model
- `s1_s2_class_comparison_summary_VV_2025.csv`: field-level Sentinel-1 summary statistics by Sentinel-2 class

Notes

- These are convenience copies of outputs generated elsewhere in the pipeline.
- Source files remain in `data/processed/validation/`.
- Raw training tables and full field-level extracts were not copied here.
