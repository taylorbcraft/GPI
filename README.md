# S2REP grassland habitat quality

This repository contains the R workflow used to classify agricultural grassland habitat quality from April 2025 Sentinel-2 imagery and field measurements in Southwest Friesland.

## Run the analysis

Open the project root in R and run:

```r
source("scripts/run_calibration_and_prediction.R")
```

This rebuilds the training data, validates S2REP against field measurements, trains the random forest, creates the field-level habitat map, compares it with Sentinel-1 LUI, and generates figures.

To apply the saved model to other years and study areas:

```r
source("scripts/run_prediction_only.R")
```

Set the prediction year and shared paths in `config.R`. Prediction rasters must follow this naming pattern:

```text
data/processed/rasters/<index>_<year>-04-median_mosaic.tif
```

## Required inputs

| Input | Purpose |
| --- | --- |
| `data/raw/environmental_field_data_2025.csv` | Soil, vegetation, and observer-assigned field measurements |
| `data/raw/plant_diversity_plots_2025.csv` | Plant richness measurements |
| `data/spatial/sampled_zone_geometry.gpkg` | Sampled training polygons |
| `data/spatial/field_geometry.gpkg` | Study-area field boundaries |
| `data/processed/rasters/*_2025-04-median_mosaic.tif` | Sentinel-2 predictor rasters |
| `data/processed/rasters/S1_VV_LogRatio_StdDev_Ascending_2025.tif` | Sentinel-1 seasonal LUI raster |
| `data/raw/biotic/` | Godwit and insect data |

The Sentinel-2 predictors are S2REP, NDVI, NDWI, SAVI, EVI, MSI, NDMI, and MNDWI.

## Workflow

| Script | Output |
| --- | --- |
| `01_build_anchor_training_data.R` | Polygon-level field and raster training data |
| `02_validate_environmental_relationships.R` | GAM results and the environmental validation figure |
| `03_define_candidate_habitat_classes.R` | KNN and rule-based candidate habitat classes |
| `04_train_habitat_classifier.R` | Cross-validated random forest and the saved model |
| `05_predict_habitat_map.R` | Field-level habitat map and map figure |
| `06_compare_s1_s2_classes.R` | Sentinel-1 comparison results and figure |
| `07_analyze_biotic_relationships.R` | Godwit and insect figures |

The optional `00_acquire_s1_imagery_rgee.R` and `00_acquire_s2_imagery_rgee.R` scripts acquire remote-sensing inputs through Google Earth Engine. They require a working `rgee` and Google Drive setup.

## Main outputs

```text
data/processed/models/habitat_best_model_2025.rds
data/processed/spatial/field_habitat_class_map_2025.gpkg
figures/environmental_validation_plots_2025.png
figures/field_habitat_class_map_2025.png
figures/s1_s2_class_comparison_VV_2025.png
figures/godwit_relationships_2025.png
figures/insect_relationships_2025.png
tables/table_1_gam_fit_statistics_2025.csv
tables/table_2_gam_terms_2025.csv
tables/table_3_random_forest_model_performance_2025.csv
tables/table_4_random_forest_class_accuracy_2025.csv
tables/table_5_random_forest_confusion_matrix_2025.csv
```

## R packages

The workflow uses `tidyverse`, `sf`, `terra`, `exactextractr`, `janitor`, `mgcv`, `caret`, `randomForest`, `patchwork`, `ggspatial`, and `ragg`. The acquisition scripts additionally require `rgee` and `googledrive`.
