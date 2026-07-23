# Apply the saved random forest and summarize pixel predictions to fields.
# The script loads the selected Sentinel-2 predictors, classifies pixels in
# memory, assigns each field its coverage-weighted dominant habitat class, and
# writes the field GeoPackage and manuscript map figure.

source("config.R")
source("scripts/shared_raster_helpers.R")

library(tidyverse)
library(randomForest)
library(sf)
library(terra)
library(exactextractr)
library(janitor)
library(ggplot2)
library(ggspatial)

# class index prediction
rf_predict_index <- function(model, data) {
  predictor_data <- as.data.frame(data)
  complete_rows <- complete.cases(predictor_data)
  predicted_class <- rep(NA_integer_, nrow(predictor_data))

  if (any(complete_rows)) {
    predicted_class[complete_rows] <- predict(
      model,
      newdata = predictor_data[complete_rows, , drop = FALSE]
    ) %>%
      as.character() %>%
      match(habitat_class_levels)
  }

  predicted_class
}

# field-level class summary
summarise_field_classes <- function(values, coverage_fraction) {
  valid_coverage <- !is.na(coverage_fraction)
  valid_values <- !is.na(values) & valid_coverage

  if (!any(valid_values)) {
    return(tibble(
      habitat_class = NA_character_,
      classified_fraction = 0,
      low_fraction = NA_real_,
      moderate_fraction = NA_real_,
      high_fraction = NA_real_
    ))
  }

  # coverage-weighted fractions
  class_fractions <- tibble(
    habitat_class = factor(
      habitat_class_levels[values[valid_values]],
      levels = habitat_class_levels
    ),
    weight = coverage_fraction[valid_values]
  ) %>%
    group_by(habitat_class) %>%
    summarise(weight = sum(weight), .groups = "drop") %>%
    complete(
      habitat_class = factor(habitat_class_levels, levels = habitat_class_levels),
      fill = list(weight = 0)
    ) %>%
    mutate(fraction = weight / sum(weight))

  tibble(
    habitat_class = class_fractions %>%
      arrange(desc(weight), habitat_class) %>%
      slice(1) %>%
      pull(habitat_class) %>%
      as.character(),
    classified_fraction = sum(class_fractions$weight) / sum(coverage_fraction[valid_coverage]),
    low_fraction = class_fractions$fraction[class_fractions$habitat_class == "low"],
    moderate_fraction = class_fractions$fraction[class_fractions$habitat_class == "moderate"],
    high_fraction = class_fractions$fraction[class_fractions$habitat_class == "high"]
  )
}

ensure_dirs(c(paths$spatial_dir, paths$figure_dir))

# load selected predictors
predictor_stack <- load_predictor_stack(
  band_names = read_csv(path_best_model_metadata(), show_col_types = FALSE)$predictors[[1]] %>%
    str_split(",\\s*") %>%
    .[[1]],
  image_date = paste0(prediction_year, "-04-median")
)

# field geometry
fields <- st_read(paths$field_geometry, quiet = TRUE) %>%
  clean_names() %>%
  rename(field_id = meadow_id) %>%
  mutate(field_id = as.character(field_id)) %>%
  st_transform(crs(predictor_stack))

# pixel prediction in memory
pixel_predictions <- terra::predict(
  object = predictor_stack,
  model = readRDS(path_best_model()),
  fun = rf_predict_index,
  overwrite = TRUE,
  wopt = list(datatype = "INT1U")
)

levels(pixel_predictions) <- data.frame(
  value = seq_along(habitat_class_levels),
  habitat_class = habitat_class_levels
)

# summarize pixel classes to fields
field_predictions <- exact_extract(
  pixel_predictions,
  fields,
  summarise_field_classes,
  progress = FALSE
) %>%
  bind_cols(fields %>% st_drop_geometry() %>% select(field_id)) %>%
  bind_cols(tibble(
    s2rep_mean = exact_extract(
      predictor_stack[["s2rep"]],
      fields,
      "mean",
      progress = FALSE
    )
  )) %>%
  relocate(field_id)

# join field predictions
field_habitat_map <- fields %>%
  left_join(field_predictions, by = "field_id")

# save field map
st_write(
  field_habitat_map,
  path_field_map(prediction_year),
  delete_dsn = TRUE,
  quiet = TRUE
)

# field map
ggsave(
  filename = file.path(paths$figure_dir, paste0("field_habitat_class_map_", prediction_year, ".png")),
  plot = ggplot(st_transform(field_habitat_map, 4326)) +
    geom_sf(aes(fill = habitat_class), color = NA) +
    scale_fill_manual(
      values = habitat_class_palette,
      breaks = habitat_class_levels,
      labels = habitat_class_levels,
      drop = FALSE,
      na.value = "grey85"
    ) +
    coord_sf(
      crs = st_crs(4326),
      datum = st_crs(4326),
      expand = FALSE
    ) +
    annotation_scale(
      location = "br",
      width_hint = 0.22,
      pad_x = grid::unit(0.35, "in"),
      pad_y = grid::unit(0.2, "in")
    ) +
    annotation_north_arrow(
      location = "tr",
      which_north = "true",
      pad_x = grid::unit(0.2, "in"),
      pad_y = grid::unit(0.2, "in"),
      style = north_arrow_orienteering
    ) +
    labs(
      x = NULL,
      y = NULL,
      fill = "Habitat quality"
    ) +
    theme_minimal(base_family = "Avenir Next", base_size = 12) +
    theme(
      text = element_text(face = "bold"),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey70", linewidth = 0.5),
      axis.text = element_text(color = "grey35"),
      legend.position = c(0.9, 0.22),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = scales::alpha("white", 0.9), color = "grey75"),
      legend.key.height = grid::unit(0.22, "in"),
      legend.title = element_text(size = 12, color = "grey20"),
      legend.text = element_text(size = 10, color = "grey35")
    ),
  device = ragg::agg_png,
  width = 8,
  height = 8,
  dpi = 300
)
