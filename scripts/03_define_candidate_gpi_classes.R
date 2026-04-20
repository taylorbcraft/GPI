# Define the supervised GPI target used for model training.
# Each sampled polygon is assigned by leave-one-out KNN in standardized field
# measurement space, using observer_estimated_GPI as the neighbor label.
#
# Input: data/processed/training/anchor_zone_training_data_<calibration_year>.csv
# Outputs: candidate training table, KNN method summary, and diagnostic boxplot.

source("config.R")

library(tidyverse)
library(patchwork)

# Field measurements used to locate each sampled polygon in ecological space.
field_feature_cols <- c(
  "soil_resistance_std",
  "soil_moisture_std",
  "plant_richness_std",
  "vegetation_height_std"
)

# Put variables on the same 0-1 scale before computing KNN distances.
min_max_rescale <- function(x) {
  rng <- range(x, na.rm = TRUE)

  if (isTRUE(all.equal(rng[1], rng[2]))) {
    return(rep(0.5, length(x)))
  }

  (x - rng[1]) / (rng[2] - rng[1])
}

# Majority voting is kept available, but the configured workflow uses weights.
majority_distance_class <- function(labels, distances) {
  counts <- table(factor(labels, levels = gpi_class_levels))
  tied_classes <- names(counts)[counts == max(counts)]

  if (length(tied_classes) == 1) {
    return(tied_classes)
  }

  mean_distances <- map_dbl(
    tied_classes,
    ~ mean(distances[labels == .x], na.rm = TRUE)
  )

  tied_classes[which.min(mean_distances)]
}

# Distance weighting lets nearer neighbors carry more influence than farther ones.
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

# Leave-one-out KNN assigns each sampled polygon from the other sampled polygons.
knn_leave_one_out <- function(data, feature_cols, k, vote_method) {
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

    if (length(candidate_rows) == 0) {
      next
    }

    distances <- sqrt(rowSums(
      sweep(feature_matrix[candidate_rows, , drop = FALSE], 2, feature_matrix[i, ], "-")^2
    ))

    nearest_order <- order(distances)[seq_len(min(k, length(candidate_rows)))]
    nearest_labels <- labels[candidate_rows[nearest_order]]
    nearest_distances <- distances[nearest_order]

    predicted[[i]] <- switch(
      vote_method,
      majority_distance = majority_distance_class(nearest_labels, nearest_distances),
      weighted_distance = weighted_distance_class(nearest_labels, nearest_distances),
      stop("Unknown knn_vote_method: ", vote_method, call. = FALSE)
    )
  }

  factor(predicted, levels = gpi_class_levels)
}

plot_target_boxplot <- function(data, target_col, title_text) {
  ggplot(data, aes(x = .data[[target_col]], y = s2rep)) +
    geom_boxplot() +
    labs(x = title_text, y = "s2rep") +
    theme_bw()
}

ensure_dirs(c(paths$training_dir, paths$validation_dir, paths$figure_dir))

require_file(path_anchor_training())

dat <- read_csv(path_anchor_training(), show_col_types = FALSE) %>%
  mutate(
    observer_estimated_GPI = factor(observer_estimated_GPI, levels = gpi_class_levels)
  )

missing_classes <- setdiff(gpi_class_levels, unique(as.character(dat$observer_estimated_GPI)))
if (length(missing_classes) > 0) {
  stop(
    paste0(
      "observer_estimated_GPI is missing required classes: ",
      paste(missing_classes, collapse = ", "),
      "."
    ),
    call. = FALSE
  )
}

# Standardize field variables, then derive the training target from KNN labels.
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
      k = knn_k,
      vote_method = knn_vote_method
    )
  )

# Keep the historical filename, but write the current KNN settings and counts.
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
      knn_vote_method,
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

# Record the target recipe beside the data so the model can be reproduced.
write_csv(knn_summary, path_estimated_thresholds())

write_csv(candidate_gpi_training_data, path_candidate_training())

# Compare observer labels with the KNN-derived target on the S2REP gradient.
p1 <- plot_target_boxplot(candidate_gpi_training_data, "observer_estimated_GPI", "observer_estimated_GPI")
p2 <- plot_target_boxplot(candidate_gpi_training_data, "gpi_class_knn_all_field", "weighted KNN all field variables")

candidate_plot <- p1 + p2 + plot_layout(ncol = 2)

ggsave(
  filename = path_candidate_class_plot(),
  plot = candidate_plot,
  width = 8,
  height = 4,
  dpi = 300
)

candidate_plot
