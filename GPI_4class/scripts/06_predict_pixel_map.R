# Apply the saved 4-class random forest to the prediction-year raster stack.
# This script writes the pixel-level map and summarizes the result
# back to the reporting polygons for field-level outputs.

source("config.R")

library(tidyverse)
library(randomForest)
library(sf)
library(terra)
library(exactextractr)
library(janitor)
library(ggplot2)

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

summarise_field_classes <- function(values, coverage_fraction) {
  coverage_valid <- !is.na(coverage_fraction)
  total_coverage <- sum(coverage_fraction[coverage_valid])
  valid <- !is.na(values) & coverage_valid
  fraction_cols <- as.list(rep(NA_real_, length(gpi_class_levels)))
  names(fraction_cols) <- paste0(gpi_class_levels, "_fraction")

  if (!any(valid)) {
    return(tibble(
      dominant_class = NA_character_,
      classified_fraction = 0,
      !!!fraction_cols
    ))
  }

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

  fraction_cols <- as.list(value_tbl$fraction)
  names(fraction_cols) <- paste0(as.character(value_tbl$gpi_class), "_fraction")

  tibble(
    dominant_class = dominant_class,
    classified_fraction = sum(value_tbl$weight) / total_coverage,
    !!!fraction_cols
  )
}

ensure_dirs(c(paths$prediction_dir, paths$spatial_dir, paths$figure_dir))

predictor_stack <- rast(path_prediction_stack())
names(predictor_stack) <- model_predictor_bands

gpi_best_model <- readRDS(path_best_model())

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

fields <- st_read(paths$field_geometry, quiet = TRUE) %>%
  clean_names() %>%
  rename(field_id = meadow_id) %>%
  mutate(field_id = as.character(field_id)) %>%
  st_transform(crs(pixel_gpi_map))

field_predictions <- exact_extract(
  pixel_gpi_map,
  fields,
  summarise_field_classes,
  progress = FALSE
) %>%
  bind_cols(fields %>% st_drop_geometry() %>% select(field_id)) %>%
  relocate(field_id)

write_csv(field_predictions, path_pixel_predictions())

field_gpi_map <- fields %>%
  left_join(field_predictions, by = "field_id")

st_write(
  field_gpi_map,
  path_field_map(),
  delete_dsn = TRUE,
  quiet = TRUE
)

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
  col = unname(gpi_class_palette[gpi_class_levels]),
  axes = FALSE,
  plg = list(x = "topright")
)

dev.off()

ggsave(
  filename = path_field_map_preview(),
  plot = ggplot(field_gpi_map) +
    geom_sf(aes(fill = dominant_class), color = NA) +
    scale_fill_manual(
      values = gpi_class_palette,
      breaks = gpi_class_levels,
      labels = unname(gpi_class_display_labels[gpi_class_levels]),
      drop = FALSE,
      na.value = "grey85"
    ) +
    theme_void() +
    labs(fill = "GPI class"),
  width = 8,
  height = 8,
  dpi = 300
)
