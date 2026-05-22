# Build the supervised class target with leave-one-out KNN.
# Observer-estimated GPI supplies the class names, and ecological similarity
# across polygons regularizes the final target used in classification.

source("config.R")

library(tidyverse)
library(patchwork)

field_feature_cols <- c(
  "soil_resistance_std",
  "soil_moisture_std",
  "plant_richness_std",
  "vegetation_height_std"
)

min_max_rescale <- function(x) {
  rng <- range(x, na.rm = TRUE)

  if (isTRUE(all.equal(rng[1], rng[2]))) {
    return(rep(0.5, length(x)))
  }

  (x - rng[1]) / (rng[2] - rng[1])
}

weighted_distance_class <- function(labels, distances) {
  weights <- 1 / (distances + 1e-6)
  weight_tbl <- tibble(label = factor(labels, levels = gpi_class_levels), weight = weights) %>%
    group_by(label) %>%
    summarise(weight = sum(weight), .groups = "drop") %>%
    complete(label = factor(gpi_class_levels, levels = gpi_class_levels), fill = list(weight = 0))

  tied_classes <- as.character(weight_tbl$label[weight_tbl$weight == max(weight_tbl$weight)])

  if (length(tied_classes) == 1) {
    return(tied_classes)
  }

  mean_distances <- map_dbl(
    tied_classes,
    ~ mean(distances[labels == .x], na.rm = TRUE)
  )

  tied_classes[which.min(mean_distances)]
}

knn_leave_one_out <- function(data, feature_cols, k) {
  feature_matrix <- data %>%
    select(all_of(feature_cols)) %>%
    as.matrix()

  labels <- as.character(data$observer_estimated_GPI)
  complete_rows <- complete.cases(feature_matrix) & !is.na(labels)
  predicted <- rep(NA_character_, nrow(data))

  for (i in seq_len(nrow(data))) {
    if (!complete_rows[[i]]) {
      next
    }

    candidate_rows <- which(complete_rows & seq_len(nrow(data)) != i)
    distances <- sqrt(rowSums(
      sweep(feature_matrix[candidate_rows, , drop = FALSE], 2, feature_matrix[i, ], "-")^2
    ))

    nearest_order <- order(distances)[seq_len(min(k, length(candidate_rows)))]
    nearest_labels <- labels[candidate_rows[nearest_order]]
    nearest_distances <- distances[nearest_order]

    predicted[[i]] <- weighted_distance_class(nearest_labels, nearest_distances)
  }

  factor(predicted, levels = gpi_class_levels)
}

ensure_dirs(c(paths$training_dir, paths$validation_dir, paths$figure_dir))

# Read anchor table and set class order
dat <- read_csv(path_anchor_training(), show_col_types = FALSE) %>%
  mutate(observer_estimated_GPI = factor(observer_estimated_GPI, levels = gpi_class_levels))

# Standardize field gradients before KNN
candidate_gpi_training_data <- dat %>%
  mutate(
    soil_resistance_std = min_max_rescale(soil_resistance),
    soil_moisture_std = min_max_rescale(soil_moisture),
    plant_richness_std = min_max_rescale(plant_richness_mean),
    vegetation_height_std = min_max_rescale(vegetation_height)
  ) %>%
  mutate(
    gpi_class_knn_all_field = knn_leave_one_out(
      data = .,
      feature_cols = field_feature_cols,
      k = knn_k
    )
  )

# Save KNN recipe and class counts
knn_summary <- bind_rows(
  tibble(
    record_type = "method",
    target = "gpi_class_knn_all_field",
    setting = c("algorithm", "k", "validation_style", "label_source", "vote_method", "tie_break"),
    value = c(
      "k_nearest_neighbors",
      as.character(knn_k),
      "leave_one_out",
      "observer_estimated_GPI",
      "weighted_distance",
      "smallest_mean_neighbor_distance"
    )
  ),
  tibble(
    record_type = "feature",
    target = "gpi_class_knn_all_field",
    setting = field_feature_cols,
    value = "min_max_standardized"
  ),
  candidate_gpi_training_data %>%
    count(gpi_class_knn_all_field, name = "n") %>%
    transmute(
      record_type = "class_count",
      target = "gpi_class_knn_all_field",
      setting = as.character(gpi_class_knn_all_field),
      value = as.character(n)
    )
)

write_csv(knn_summary, path_estimated_thresholds())
write_csv(candidate_gpi_training_data, path_candidate_training())

# Compare observer and KNN labels on S2REP
p1 <- ggplot(candidate_gpi_training_data, aes(x = observer_estimated_GPI, y = s2rep)) +
  geom_boxplot() +
  labs(x = "observer_estimated_GPI", y = "s2rep") +
  theme_bw()

p2 <- ggplot(candidate_gpi_training_data, aes(x = gpi_class_knn_all_field, y = s2rep)) +
  geom_boxplot() +
  labs(x = "weighted KNN all field variables", y = "s2rep") +
  theme_bw()

candidate_plot <- p1 + p2 + plot_layout(ncol = 2)

ggsave(
  filename = path_candidate_class_plot(),
  plot = candidate_plot,
  width = 8,
  height = 4,
  dpi = 300
)

candidate_plot
