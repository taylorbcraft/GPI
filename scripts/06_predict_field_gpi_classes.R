# Apply the saved GPI model to the full field layer.
# Produces the prediction table, joins predictions back to field geometry, and
# writes the final geopackage plus a preview map.
#
# Inputs: field predictor table, saved model RDS, and field_geometry.gpkg.
# Outputs: field predictions CSV, field GPI geopackage, and map preview PNG.

source("config.R")

library(tidyverse)
library(randomForest)
library(sf)
library(janitor)
library(ggplot2)

ensure_dirs(c(paths$prediction_dir, paths$spatial_dir, paths$figure_dir))

require_file(path_field_predictor_data())
require_file(path_best_model())
require_file(paths$field_geometry)

field_predictor_data <- read_csv(path_field_predictor_data(), show_col_types = FALSE)

gpi_best_model <- readRDS(path_best_model())

# Predict only fields with complete values for the model predictor set.
field_prediction_data <- field_predictor_data %>%
  select(field_id, all_of(model_predictor_bands)) %>%
  drop_na(all_of(model_predictor_bands))

field_predictions <- field_prediction_data %>%
  mutate(
    gpi_class = predict(gpi_best_model, newdata = select(., all_of(model_predictor_bands))),
    gpi_class = factor(gpi_class, levels = gpi_class_levels)
  )

write_csv(field_predictions, path_field_predictions())

# Join predictions back to geometry so the tabular output and map stay aligned.
fields <- st_read(paths$field_geometry, quiet = TRUE) %>%
  clean_names()

if (!field_id_col %in% names(fields)) {
  stop(paste0("Field identifier column not found in field geometry: ", field_id_col), call. = FALSE)
}

names(fields)[names(fields) == field_id_col] <- "field_id"

fields <- fields %>%
  mutate(field_id = as.character(field_id))

field_gpi_map <- fields %>%
  left_join(
    field_predictions %>%
      mutate(field_id = as.character(field_id)),
    by = "field_id"
  )

st_write(
  field_gpi_map,
  path_field_map(),
  delete_dsn = TRUE,
  quiet = TRUE
)

# Fixed palette keeps class colors stable across annual maps.
gpi_map_palette <- c(
  extensive = "#2E7D32",
  mid = "#FDAE61",
  intensive = "#D73027"
)

p <- ggplot(field_gpi_map) +
  geom_sf(aes(fill = gpi_class), color = NA) +
  scale_fill_manual(
    values = gpi_map_palette,
    breaks = gpi_class_levels,
    labels = c("Extensive", "Mid", "Intensive"),
    drop = FALSE,
    na.value = "grey85"
  ) +
  theme_void() +
  labs(fill = "GPI class")

ggsave(
  filename = path_field_map_preview(),
  plot = p,
  width = 8,
  height = 8,
  dpi = 300
)

p
