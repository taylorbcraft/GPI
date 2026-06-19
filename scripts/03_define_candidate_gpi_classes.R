# Build supervised class targets with KNN and RBC

source("config.R")

library(tidyverse)
library(patchwork)
library(caret)

knn_k <- 7

field_feature_cols <- c(
  "soil_resistance_std",
  "soil_moisture_std",
  "plant_richness_std",
  "vegetation_height_std"
)

# Put all field variables on a common  scale before distance-based
# classification so no single variable dominates because of its units
min_max_rescale <- function(x) {
  rng <- range(x, na.rm = TRUE)

  if (isTRUE(all.equal(rng[1], rng[2]))) {
    return(rep(0.5, length(x)))
  }

  (x - rng[1]) / (rng[2] - rng[1])
}

# Summarize how closely each derived target reproduces the observer labels
classification_metrics <- function(reference, prediction) {
  cm <- confusionMatrix(
    data = factor(prediction, levels = gpi_class_levels),
    reference = factor(reference, levels = gpi_class_levels)
  )

  tibble(
    accuracy = unname(cm$overall["Accuracy"]),
    kappa = unname(cm$overall["Kappa"])
  )
}

# Convert the neighbor set into one class label using inverse-distance
# weighting. If two classes tie on total weight, prefer the one with the
# smaller mean neighbor distance
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

# Leave-one-out KNN assigns each polygon using all other polygons as the
# reference set, which avoids letting a polygon vote for its own label
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
    # Distances are measured in standardized field-variable space
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

# Build a rule-based comparison target from the standardized soil
# resistance and richness gradients. The cutoffs are anchored to the
# observer-defined extensive and intensive class distributions
derive_rbc_3class <- function(data) {
  valid <- complete.cases(data[, c("soil_resistance_std", "plant_richness_std", "observer_estimated_GPI")])
  work <- data[valid, ]

  class_iqr <- work %>%
    group_by(observer_estimated_GPI) %>%
    summarise(
      resistance_q1 = quantile(soil_resistance_std, probs = 0.25, na.rm = TRUE),
      resistance_q3 = quantile(soil_resistance_std, probs = 0.75, na.rm = TRUE),
      richness_q1 = quantile(plant_richness_std, probs = 0.25, na.rm = TRUE),
      richness_q3 = quantile(plant_richness_std, probs = 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(observer_estimated_GPI = factor(observer_estimated_GPI, levels = gpi_class_levels)) %>%
    arrange(observer_estimated_GPI)

  extensive_iqr <- class_iqr %>% filter(observer_estimated_GPI == "extensive")
  intensive_iqr <- class_iqr %>% filter(observer_estimated_GPI == "intensive")

  resistance_cutoff <- mean(c(extensive_iqr$resistance_q3, intensive_iqr$resistance_q1))
  richness_cutoff <- mean(c(extensive_iqr$richness_q1, intensive_iqr$richness_q3))

  # Assign the extremes when both variables support the same direction;
  # otherwise route ambiguous cases to the middle class
  predicted_all <- rep(NA_character_, nrow(data))
  predicted_all[valid] <- case_when(
    work$soil_resistance_std <= resistance_cutoff & work$plant_richness_std >= richness_cutoff ~ "extensive",
    work$soil_resistance_std > resistance_cutoff & work$plant_richness_std < richness_cutoff ~ "intensive",
    TRUE ~ "mid"
  )

  metrics <- classification_metrics(work$observer_estimated_GPI, predicted_all[valid])

  list(
    classes = factor(predicted_all, levels = gpi_class_levels),
    summary = bind_rows(
      tibble(
        record_type = "method",
        target = "gpi_class_rbc_soil_richness",
        setting = c("algorithm", "resistance_cutoff", "richness_cutoff", "input_variables", "threshold_source", "threshold_scale", "label_reference", "agreement_accuracy", "agreement_kappa"),
        value = c(
          "rule_based_classification",
          format(resistance_cutoff, digits = 6),
          format(richness_cutoff, digits = 6),
          "soil_resistance_std, plant_richness_std",
          "initial_observer_class_iqr_ranges",
          "min_max_standardized",
          "observer_estimated_GPI",
          format(metrics$accuracy, digits = 6),
          format(metrics$kappa, digits = 6)
        )
      ),
      class_iqr %>%
        mutate(
          target = "gpi_class_rbc_soil_richness",
          record_type = "class_iqr",
          class = as.character(observer_estimated_GPI)
        ) %>%
        pivot_longer(
          cols = c(resistance_q1, resistance_q3, richness_q1, richness_q3),
          names_to = "setting",
          values_to = "value"
        ) %>%
        transmute(
          record_type,
          target,
          setting = paste(class, setting, sep = "_"),
          value = format(value, digits = 6)
        )
    )
  )
}

ensure_dirs(c(paths$training_dir, paths$validation_dir, paths$figure_dir))

# Read anchor table and set class order
dat <- read_csv(path_anchor_training(), show_col_types = FALSE) %>%
  mutate(observer_estimated_GPI = factor(observer_estimated_GPI, levels = gpi_class_levels))

# Standardize the field gradients once, then use the same scaled variables
# for both the KNN target and the rule-based comparison target.
candidate_gpi_training_data <- dat %>%
  mutate(
    soil_resistance_std = min_max_rescale(soil_resistance),
    soil_moisture_std = min_max_rescale(soil_moisture),
    plant_richness_std = min_max_rescale(plant_richness_mean),
    vegetation_height_std = min_max_rescale(vegetation_height)
  )

candidate_gpi_training_data <- candidate_gpi_training_data %>%
  mutate(
    gpi_class_knn_all_field = knn_leave_one_out(
      data = .,
      feature_cols = field_feature_cols,
      k = knn_k
    )
  )

rbc_result <- derive_rbc_3class(candidate_gpi_training_data)

candidate_gpi_training_data <- candidate_gpi_training_data %>%
  mutate(gpi_class_rbc_soil_richness = rbc_result$classes)

# Write table that records the KNN setup, the RBC thresholds, and
# the resulting class counts for both targets
target_summary <- bind_rows(
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
    ),
  rbc_result$summary,
  candidate_gpi_training_data %>%
    count(gpi_class_rbc_soil_richness, name = "n") %>%
    transmute(
      record_type = "class_count",
      target = "gpi_class_rbc_soil_richness",
      setting = as.character(gpi_class_rbc_soil_richness),
      value = as.character(n)
    )
)

write_csv(target_summary, path_estimated_thresholds())
write_csv(candidate_gpi_training_data, path_candidate_training())

# A quick visual check of whether each target produces the expected ordering
# along the continuous S2REP gradient
p1 <- ggplot(candidate_gpi_training_data, aes(x = observer_estimated_GPI, y = s2rep)) +
  geom_boxplot() +
  labs(x = "observer_estimated_GPI", y = "s2rep") +
  theme_bw()

p2 <- ggplot(candidate_gpi_training_data, aes(x = gpi_class_knn_all_field, y = s2rep)) +
  geom_boxplot() +
  labs(x = "weighted KNN", y = "s2rep") +
  theme_bw()

p3 <- ggplot(candidate_gpi_training_data, aes(x = gpi_class_rbc_soil_richness, y = s2rep)) +
  geom_boxplot() +
  labs(x = "RBC soil + richness", y = "s2rep") +
  theme_bw()

candidate_plot <- p1 + p2 + p3 + plot_layout(ncol = 3)

ggsave(
  filename = path_candidate_class_plot(),
  plot = candidate_plot,
  width = 12,
  height = 4,
  dpi = 300
)

candidate_plot
