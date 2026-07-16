# Build the prediction-year raster stack and apply the saved random forest.
# Write the pixel-level GPI map, summarize predicted classes to fields,
# and save preview figures for quick checks.

source("config.R")

library(tidyverse)
library(randomForest)
library(sf)
library(terra)
library(exactextractr)
library(janitor)
library(ggplot2)

# class index prediction
rf_predict_index <- function(model, data) {
  predictor_df <- as.data.frame(data)
  complete_rows <- complete.cases(predictor_df)
  out <- rep(NA_integer_, nrow(predictor_df))

  if (any(complete_rows)) {
    pred <- predict(model, newdata = predictor_df[complete_rows, , drop = FALSE])
    out[complete_rows] <- match(as.character(pred), gpi_class_levels)
  }

  out
}

# field-level class summary
summarise_field_classes <- function(values, coverage_fraction) {
  coverage_valid <- !is.na(coverage_fraction)
  total_coverage <- sum(coverage_fraction[coverage_valid])
  valid <- !is.na(values) & coverage_valid

  if (!any(valid)) {
    return(tibble(
      dominant_class = NA_character_,
      classified_fraction = 0,
      extensive_fraction = NA_real_,
      mid_fraction = NA_real_,
      intensive_fraction = NA_real_
    ))
  }

  # coverage-weighted class fractions
  value_tbl <- tibble(
    class_id = values[valid],
    weight = coverage_fraction[valid]
  ) %>%
    mutate(gpi_class = factor(gpi_class_levels[class_id], levels = gpi_class_levels)) %>%
    group_by(gpi_class) %>%
    summarise(weight = sum(weight), .groups = "drop") %>%
    complete(gpi_class = factor(gpi_class_levels, levels = gpi_class_levels), fill = list(weight = 0)) %>%
    mutate(fraction = weight / sum(weight))

  dominant_class <- value_tbl %>%
    arrange(desc(weight), gpi_class) %>%
    slice(1) %>%
    pull(gpi_class) %>%
    as.character()

  tibble(
    dominant_class = dominant_class,
    classified_fraction = sum(value_tbl$weight) / total_coverage,
    extensive_fraction = value_tbl$fraction[value_tbl$gpi_class == "extensive"],
    mid_fraction = value_tbl$fraction[value_tbl$gpi_class == "mid"],
    intensive_fraction = value_tbl$fraction[value_tbl$gpi_class == "intensive"]
  )
}

ensure_dirs(c(paths$prediction_dir, paths$spatial_dir, paths$figure_dir))

# load saved model
gpi_best_model <- readRDS(path_best_model())
best_model_metadata <- read_csv(path_best_model_metadata(), show_col_types = FALSE)
training_unit <- best_model_metadata$training_unit[[1]]
selected_predictors <- best_model_metadata$predictors[[1]] %>%
  str_split(",\\s*") %>%
  .[[1]]

if (!training_unit %in% c("pixel", "field")) {
  stop("Saved model training_unit must be 'pixel' or 'field'.")
}

if (length(selected_predictors) == 0 || any(is.na(selected_predictors))) {
  stop("Saved model metadata must include at least one predictor band.")
}

# build prediction stack
predictor_stack <- load_predictor_stack(
  band_names = selected_predictors,
  image_date = prediction_image_date
)

writeRaster(
  predictor_stack,
  filename = path_prediction_stack(),
  overwrite = TRUE
)

fields <- st_read(paths$field_geometry, quiet = TRUE) %>%
  clean_names() %>%
  rename(field_id = meadow_id) %>%
  mutate(field_id = as.character(field_id)) %>%
  st_transform(crs(predictor_stack))

if (training_unit == "pixel") {
  # predict pixel classes
  pixel_gpi_map <- terra::predict(
    object = predictor_stack,
    model = gpi_best_model,
    fun = rf_predict_index,
    overwrite = TRUE,
    wopt = list(datatype = "INT1U")
  )

  if (file.exists(path_pixel_map())) {
    file.remove(path_pixel_map())
  }

  writeRaster(
    pixel_gpi_map,
    filename = path_pixel_map(),
    overwrite = TRUE,
    datatype = "INT1U",
    gdal = c("BIGTIFF=YES", "COMPRESS=LZW")
  )

  levels(pixel_gpi_map) <- data.frame(
    value = seq_along(gpi_class_levels),
    gpi_class = gpi_class_levels
  )

  # summarize fields from raster
  field_predictions <- exact_extract(
    pixel_gpi_map,
    fields,
    summarise_field_classes,
    progress = FALSE
  ) %>%
    bind_cols(fields %>% st_drop_geometry() %>% select(field_id)) %>%
    relocate(field_id)

  field_gpi_map <- fields %>%
    left_join(field_predictions, by = "field_id")
} else {
  # summarize predictors to fields
  field_predictors <- exact_extract(
    predictor_stack,
    fields,
    "mean",
    progress = FALSE
  ) %>%
    as_tibble() %>%
    rename_with(~ str_remove(.x, "^mean\\."))

  predictor_complete <- complete.cases(field_predictors[, selected_predictors, drop = FALSE])
  predicted_class <- rep(NA_character_, nrow(field_predictors))

  if (any(predictor_complete)) {
    predicted_class[predictor_complete] <- predict(
      gpi_best_model,
      newdata = field_predictors[predictor_complete, selected_predictors, drop = FALSE]
    ) %>%
      as.character()
  }

  field_predictions <- fields %>%
    st_drop_geometry() %>%
    select(field_id) %>%
    bind_cols(field_predictors) %>%
    mutate(
      dominant_class = predicted_class,
      classified_fraction = if_else(predictor_complete, 1, 0),
      extensive_fraction = if_else(dominant_class == "extensive", 1, 0, missing = NA_real_),
      mid_fraction = if_else(dominant_class == "mid", 1, 0, missing = NA_real_),
      intensive_fraction = if_else(dominant_class == "intensive", 1, 0, missing = NA_real_)
    )

  field_gpi_map <- fields %>%
    left_join(field_predictions, by = "field_id")

  pixel_gpi_map <- rasterize(
    x = vect(field_gpi_map),
    y = predictor_stack[[1]],
    field = match(field_gpi_map$dominant_class, gpi_class_levels),
    background = NA
  )

  writeRaster(
    pixel_gpi_map,
    filename = path_pixel_map(),
    overwrite = TRUE,
    datatype = "INT1U",
    gdal = c("BIGTIFF=YES", "COMPRESS=LZW")
  )

  levels(pixel_gpi_map) <- data.frame(
    value = seq_along(gpi_class_levels),
    gpi_class = gpi_class_levels
  )
}

write_csv(field_predictions, path_pixel_predictions())

st_write(
  field_gpi_map,
  path_field_map(),
  delete_dsn = TRUE,
  quiet = TRUE
)

# save preview maps
preview_raster <- aggregate(pixel_gpi_map, fact = 4, fun = modal, na.rm = TRUE)

png(
  filename = path_pixel_map_preview(),
  width = 2400,
  height = 2400,
  res = 300
)

plot(
  preview_raster,
  type = "classes",
  col = c("#2E7D32", "#FDAE61", "#D73027"),
  axes = FALSE,
  plg = list(x = "topright")
)

dev.off()

# save field preview
ggsave(
  filename = path_field_map_preview(),
  plot = ggplot(field_gpi_map) +
    geom_sf(aes(fill = dominant_class), color = NA) +
    scale_fill_manual(
      values = c(extensive = "#2E7D32", mid = "#FDAE61", intensive = "#D73027"),
      breaks = gpi_class_levels,
      labels = c("Extensive", "Mid", "Intensive"),
      drop = FALSE,
      na.value = "grey85"
    ) +
    theme_void() +
    labs(fill = "GPI class"),
  width = 8,
  height = 8,
  dpi = 300
)
