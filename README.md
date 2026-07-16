# GPI Workflow

This project builds a reproducible workflow for mapping grassland production intensity (GPI)
from Sentinel-2 imagery and field-calibrated environmental measures. The model is calibrated
with 2025 field surveys, where sampled zones were assigned observer-based management classes and
measured for soil resistance, soil moisture, vegetation height, and plant richness. These field
measurements are used to derive a three-class training dataset using KNN. A random forest
classifier is then trained on pixels extracted from the labeled habitat polygons and applied to
the full raster stack to produce a pixel-level GPI map.

Once calibrated, the saved model can be applied to other years without collecting new field data,
provided that matching Sentinel-2 predictor rasters are available for the target year.

## Quick Start

1. Check `config.R` for the calibration year, prediction year, image dates,
   class order, and expected predictor bands.
2. To rebuild the calibration model and make the configured prediction map, run:

```r
source("scripts/run_calibration_and_prediction.R")
```

3. To apply the existing calibration model to a new prediction year, update
   `prediction_year` and `prediction_image_date`, add those rasters, then run:

```r
source("scripts/run_prediction_only.R")
```

## Workflow Logic

Sampled zones carry field measurements and observer labels. Those labels are converted into a
weighted KNN-derived GPI target. The labeled zones are then used as training masks, so all
predictor pixels inside each zone inherit the zone label. The trained random forest is applied
pixel by pixel to the prediction-year raster stack, and the resulting raster is summarized back
to fields for polygon-based reporting.

| Step | Script | Role | Main outputs |
| --- | --- | --- | --- |
| 00 | `scripts/00_export_remote_sensing_from_gee.js` | Documents the raster export used by the R workflow. Run only when rasters need to be rebuilt. | `<band>_<date>_mosaic.tif` |
| 01 | `scripts/01_build_anchor_training_data.R` | Joins sampled-zone geometry, field observations, plant richness, and calibration raster summaries into the anchor training table. | `anchor_zone_training_data_<calibration_year>.csv` |
| 02 | `scripts/02_validate_environmental_relationships.R` | Checks whether `s2rep` has interpretable relationships with measured ecological variables. | validation summary CSV, coefficients CSV, full model summary text, and figure |
| 03 | `scripts/03_define_candidate_gpi_classes.R` | Defines the supervised target with leave-one-out KNN over standardized field variables. | `candidate_gpi_training_data_<calibration_year>.csv`, KNN summary, boxplot |
| 04 | `scripts/04_train_gpi_classifier.R` | Extracts labeled pixels from sampled zones, runs forward predictor selection with grouped cross-validation, and saves diagnostics plus the selected model. | model RDS and validation tables |
| 05 | `scripts/05_predict_pixel_map.R` | Builds the prediction-year predictor raster stack, applies the saved calibration model pixel by pixel, summarizes the classified raster back to fields, and writes the final map products. | `pixel_predictor_stack_<prediction_year>.tif`, `pixel_gpi_map_<prediction_year>.tif`, `field_gpi_pixel_summary_<prediction_year>.csv`, `field_gpi_map_<prediction_year>.gpkg`, PNG previews |

## Project Layout

```text
GPI_Project/
├── config.R
├── README.md
├── data/
│   ├── raw/                 # sampled field inputs
│   ├── spatial/             # sampled-zone and full-field geometry
│   └── processed/
│       ├── models/          # fitted random forest model
│       ├── predictions/     # prediction stacks and field summaries
│       ├── rasters/         # Sentinel-2 predictor stack
│       ├── spatial/         # pixel maps and field geopackages
│       ├── training/        # anchor and candidate training tables
│       └── validation/      # diagnostics and model evaluation tables
├── figures/                 # validation and map preview figures
└── scripts/
```

## Inputs

`data/raw/environmental_field_data_<calibration_year>.csv`

Contains sampled-zone field measurements such as soil moisture, soil resistance,
vegetation height, and the observer label.

`data/raw/plant_diversity_plots_<calibration_year>.csv`

Contains plot-level plant richness observations that are summarized to `polygon_id`.

`data/spatial/sampled_zone_geometry.gpkg`

Contains the sampled training polygons used for raster extraction and joins to field observations.

`data/spatial/field_geometry.gpkg`

Contains the geometry of all fields in the study area.

`data/processed/rasters/`

Stores the predictor stack expected by calibration and prediction:

- `s2rep_<date>_mosaic.tif`
- `ndvi_<date>_mosaic.tif`
- `ndwi_<date>_mosaic.tif`
- `savi_<date>_mosaic.tif`
- `evi_<date>_mosaic.tif`
- `msi_<date>_mosaic.tif`
- `ndmi_<date>_mosaic.tif`
- `mndwi_<date>_mosaic.tif`

The calibration rasters use `calibration_image_date`. Annual maps use
`prediction_image_date`. Keep the band names and index formulas consistent
between calibration and prediction years.

## Outputs

Training outputs:

- `data/processed/training/anchor_zone_training_data_<calibration_year>.csv`
- `data/processed/training/candidate_gpi_training_data_<calibration_year>.csv`

Validation outputs:

- environmental validation summary, coefficients, full model summaries, and figure
- KNN method summary and class counts
- random forest confusion matrix, class accuracy, variable importance,
  model comparison, and selected-model metadata

Model output:

- `data/processed/models/gpi_best_model_<calibration_year>.rds`

Prediction and map outputs:

- `data/processed/predictions/pixel_predictor_stack_<prediction_year>.tif`
- `data/processed/predictions/field_gpi_pixel_summary_<prediction_year>.csv`
- `data/processed/spatial/pixel_gpi_map_<prediction_year>.tif`
- `data/processed/spatial/field_gpi_map_<prediction_year>.gpkg`
- `figures/pixel_gpi_map_<prediction_year>.png`
- `figures/field_gpi_map_<prediction_year>.png`

## Configuration

`config.R` centralizes the values that must stay consistent across scripts:

- `calibration_year` and `calibration_image_date`
- `prediction_year` and `prediction_image_date`
- expected predictor bands
- candidate predictor bands and maximum predictors retained
- ordered three-class GPI levels
- canonical input and output paths

Update `prediction_year` and `prediction_image_date` for annual mapping. Update
`predictor_bands` only when the raster stack changes.

## Required R Packages

- `tidyverse`
- `sf`
- `terra`
- `exactextractr`
- `janitor`
- `mgcv`
- `broom`
- `patchwork`
- `caret`
- `randomForest`
