# GPI Workflow

This project calibrates a GPI classifier from the 2025 field campaign and applies
that saved model to field-level Sentinel-2 summaries. The local R workflow starts
after the predictor rasters have been exported; the Earth Engine script is kept as
the record for rebuilding those rasters when needed.

## Quick Start

1. Check `config.R` for the calibration year, prediction year, image dates,
   field id column, class order, and expected predictor bands.
2. Confirm the required inputs are present:
   - `data/raw/environmental_field_data_<calibration_year>.csv`
   - `data/raw/plant_diversity_plots_<calibration_year>.csv`
   - `data/spatial/sampled_zone_geometry.gpkg`
   - `data/spatial/field_geometry.gpkg`
   - calibration rasters named `<band>_<calibration_image_date>_mosaic.tif`
   - prediction rasters named `<band>_<prediction_image_date>_mosaic.tif`
3. To rebuild the calibration model and make the configured prediction map, run:

```r
source("scripts/run_calibration_and_prediction.R")
```

4. To apply the existing calibration model to a new prediction year, update
   `prediction_year` and `prediction_image_date`, add those rasters, then run:

```r
source("scripts/run_prediction_only.R")
```

The full runner sources scripts `01` through `06`. The prediction runner sources
only scripts `05` and `06`, so it does not require new field observations.

## Workflow Logic

The project uses two spatial units:

- `zone`: sampled training polygons identified by `polygon_id`
- `field`: the full mapping layer that receives final predictions

Sampled zones carry field measurements and observer labels. The observer labels
are collapsed to three classes: original `extensive` and `mid_low` become
`extensive`, original `mid_high` becomes `mid`, and original `intensive` remains
`intensive`. Those labels are converted into a weighted KNN-derived GPI target
and used to train a random forest. The trained model is then applied to every
mapped field using the same raster summaries.

| Step | Script | Role | Main outputs |
| --- | --- | --- | --- |
| 00 | `scripts/00_export_remote_sensing_from_gee.js` | Documents the Sentinel-2 raster export used by the R workflow. Run only when rasters need to be rebuilt. | `<band>_<date>_mosaic.tif` |
| 01 | `scripts/01_build_anchor_training_data.R` | Joins sampled-zone geometry, 2025 field observations, plant richness, and calibration raster summaries into the anchor training table. | `anchor_zone_training_data_<calibration_year>.csv` |
| 02 | `scripts/02_validate_environmental_relationships.R` | Checks whether `s2rep` has interpretable relationships with measured ecological variables. | validation summary CSV and plot |
| 03 | `scripts/03_define_candidate_gpi_classes.R` | Defines the supervised target with leave-one-out KNN over standardized field variables. | `candidate_gpi_training_data_<calibration_year>.csv`, KNN summary, boxplot |
| 04 | `scripts/04_train_gpi_classifier.R` | Tunes and fits the random forest classifier, then saves model diagnostics and the selected model. | model RDS and validation tables |
| 05 | `scripts/05_build_field_prediction_data.R` | Extracts prediction-year model raster summaries for every field polygon. | `field_predictor_data_<prediction_year>.csv` |
| 06 | `scripts/06_predict_field_gpi_classes.R` | Applies the saved calibration model to field predictors and writes the final table, geopackage, and preview map. | `field_gpi_predictions_<prediction_year>.csv`, `field_gpi_map_<prediction_year>.gpkg`, PNG preview |

## Project Layout

```text
GPI_Project/
├── config.R
├── README.md
├── data/
│   ├── raw/                 # non-spatial sampled field inputs
│   ├── spatial/             # sampled-zone and full-field geometry
│   └── processed/
│       ├── models/          # fitted random forest model
│       ├── predictions/     # field-level prediction tables
│       ├── rasters/         # Sentinel-2 predictor stack
│       ├── spatial/         # final mapped GPI output
│       ├── training/        # anchor and candidate training tables
│       └── validation/      # diagnostics and model evaluation tables
├── figures/                 # validation and map preview figures
└── scripts/
```

## Inputs

`data/raw/environmental_field_data_<calibration_year>.csv`

Contains sampled-zone field measurements such as soil moisture, soil resistance,
vegetation height, and the observer label. The source column `in_lui` is renamed
to `observer_estimated_GPI` in script `01`.

`data/raw/plant_diversity_plots_<calibration_year>.csv`

Contains plot-level plant richness observations that are summarized to
`polygon_id`.

`data/spatial/sampled_zone_geometry.gpkg`

Defines the sampled training polygons used for raster extraction and joins to
field observations.

`data/spatial/field_geometry.gpkg`

Defines the full mapping layer. The configured `field_id_col` is standardized to
`field_id` while the prediction tables are being built.

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

- environmental validation summary and plot
- KNN method summary and class counts
- random forest tuning, confusion matrix, class accuracy, variable importance,
  model comparison, and selected-model metadata

Model output:

- `data/processed/models/gpi_best_model_<calibration_year>.rds`

Prediction and map outputs:

- `data/processed/predictions/field_predictor_data_<prediction_year>.csv`
- `data/processed/predictions/field_gpi_predictions_<prediction_year>.csv`
- `data/processed/spatial/field_gpi_map_<prediction_year>.gpkg`
- `figures/field_gpi_map_<prediction_year>.png`

Note: `gpi_estimated_rule_thresholds_<calibration_year>.csv` is a compatibility
filename. In the current workflow it stores KNN settings and class counts, not
numeric rule thresholds.

## Configuration

`config.R` centralizes the values that must stay consistent across scripts:

- `calibration_year` and `calibration_image_date`
- `prediction_year` and `prediction_image_date`
- geometry id columns
- expected predictor bands
- model predictor bands
- ordered three-class GPI levels
- KNN settings
- canonical input and output paths

Update `prediction_year` and `prediction_image_date` for annual mapping. Update
`predictor_bands` only when the raster stack changes. Update `field_id_col` only
when the full field geometry uses a different source id column.

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
