# Build candidate supervised habitat-class targets from the field data.
# The script rescales the field measurements, creates leave-one-out weighted KNN
# labels, derives a rule-based soil-resistance and plant-richness target, and
# compares S2REP ordering across the candidate classifications.

source("config.R")

library(tidyverse)
library(patchwork)

# rescale field gradients
min_max_rescale <- function(x) {
  rng <- range(x, na.rm = TRUE)

  if (isTRUE(all.equal(rng[1], rng[2]))) {
    return(rep(0.5, length(x)))
  }

  (x - rng[1]) / (rng[2] - rng[1])
}

# inverse-distance vote
weighted_distance_class <- function(labels, distances) {
  weights <- 1 / (distances + 1e-6)
  weight_tbl <- tibble(label = factor(labels, levels = habitat_class_levels), weight = weights) %>%
    group_by(label) %>%
    summarise(weight = sum(weight), .groups = "drop") %>%
    complete(label = factor(habitat_class_levels, levels = habitat_class_levels), fill = list(weight = 0))

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

# leave-one-out knn labels
knn_leave_one_out <- function(data, feature_cols, k) {
  feature_matrix <- data %>%
    select(all_of(feature_cols)) %>%
    as.matrix()

  labels <- as.character(data$observer_estimated_class)
  complete_rows <- complete.cases(feature_matrix) & !is.na(labels)
  predicted <- rep(NA_character_, nrow(data))

  for (i in seq_len(nrow(data))) {
    if (!complete_rows[[i]]) {
      next
    }

    candidate_rows <- which(complete_rows & seq_len(nrow(data)) != i)
    # distances in scaled field space
    distances <- sqrt(rowSums(
      sweep(feature_matrix[candidate_rows, , drop = FALSE], 2, feature_matrix[i, ], "-")^2
    ))

    nearest_order <- order(distances)[seq_len(min(k, length(candidate_rows)))]
    nearest_labels <- labels[candidate_rows[nearest_order]]
    nearest_distances <- distances[nearest_order]

    predicted[[i]] <- weighted_distance_class(nearest_labels, nearest_distances)
  }

  factor(predicted, levels = habitat_class_levels)
}

# rule-based comparison target
derive_rbc_3class <- function(data) {
  valid <- complete.cases(data[, c("soil_resistance_std", "plant_richness_std", "observer_estimated_class")])
  work <- data[valid, ]

  class_iqr <- work %>%
    group_by(observer_estimated_class) %>%
    summarise(
      resistance_q1 = quantile(soil_resistance_std, probs = 0.25, na.rm = TRUE),
      resistance_q3 = quantile(soil_resistance_std, probs = 0.75, na.rm = TRUE),
      richness_q1 = quantile(plant_richness_std, probs = 0.25, na.rm = TRUE),
      richness_q3 = quantile(plant_richness_std, probs = 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(observer_estimated_class = factor(observer_estimated_class, levels = habitat_class_levels)) %>%
    arrange(observer_estimated_class)

  high_iqr <- class_iqr %>% filter(observer_estimated_class == "high")
  low_iqr <- class_iqr %>% filter(observer_estimated_class == "low")

  resistance_cutoff <- mean(c(high_iqr$resistance_q3, low_iqr$resistance_q1))
  richness_cutoff <- mean(c(high_iqr$richness_q1, low_iqr$richness_q3))

  # threshold-based class assignment
  predicted_all <- rep(NA_character_, nrow(data))
  predicted_all[valid] <- case_when(
    work$soil_resistance_std <= resistance_cutoff & work$plant_richness_std >= richness_cutoff ~ "high",
    work$soil_resistance_std > resistance_cutoff & work$plant_richness_std < richness_cutoff ~ "low",
    TRUE ~ "moderate"
  )

  factor(predicted_all, levels = habitat_class_levels)
}

ensure_dirs(c(paths$training_dir, paths$figure_dir))

# anchor training data
dat <- read_csv(path_anchor_training(), show_col_types = FALSE) %>%
  mutate(observer_estimated_class = factor(observer_estimated_class, levels = habitat_class_levels))

# standardized field variables
candidate_habitat_training_data <- dat %>%
  mutate(
    soil_resistance_std = min_max_rescale(soil_resistance),
    soil_moisture_std = min_max_rescale(soil_moisture),
    plant_richness_std = min_max_rescale(plant_richness_mean),
    vegetation_height_std = min_max_rescale(vegetation_height)
  )

candidate_habitat_training_data <- candidate_habitat_training_data %>%
  mutate(
    habitat_class_knn_all_field = knn_leave_one_out(
      data = .,
      feature_cols = c(
        "soil_resistance_std",
        "soil_moisture_std",
        "plant_richness_std",
        "vegetation_height_std"
      ),
      k = 7
    )
  )

candidate_habitat_training_data <- candidate_habitat_training_data %>%
  mutate(habitat_class_rbc_soil_richness = derive_rbc_3class(.))

# save candidate labels
write_csv(candidate_habitat_training_data, path_candidate_training())

# class ordering boxplots
p1 <- ggplot(
  candidate_habitat_training_data,
  aes(x = observer_estimated_class, y = s2rep, fill = observer_estimated_class)
) +
  geom_boxplot(linewidth = 0.8, color = "grey20", alpha = 0.9) +
  scale_fill_manual(values = habitat_class_palette, drop = FALSE) +
  labs(x = "Habitat quality", y = "S2REP", title = "Observer estimate")

p2 <- ggplot(
  candidate_habitat_training_data,
  aes(x = habitat_class_knn_all_field, y = s2rep, fill = habitat_class_knn_all_field)
) +
  geom_boxplot(linewidth = 0.8, color = "grey20", alpha = 0.9) +
  scale_fill_manual(values = habitat_class_palette, drop = FALSE) +
  labs(x = "Habitat quality", y = "S2REP", title = "Weighted KNN")

p3 <- ggplot(
  candidate_habitat_training_data,
  aes(x = habitat_class_rbc_soil_richness, y = s2rep, fill = habitat_class_rbc_soil_richness)
) +
  geom_boxplot(linewidth = 0.8, color = "grey20", alpha = 0.9) +
  scale_fill_manual(values = habitat_class_palette, drop = FALSE) +
  labs(x = "Habitat quality", y = "S2REP", title = "Rule-based classes")

candidate_plot <- (
  p1 + p2 + p3 +
    plot_layout(ncol = 3, axes = "collect_y", axis_titles = "collect")
) &
  theme_minimal(base_family = "Avenir Next", base_size = 18) &
  theme(
    text = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.4),
    axis.text = element_text(color = "grey35"),
    axis.title = element_text(color = "grey20"),
    plot.title = element_text(face = "bold", size = 18, color = "grey20"),
    plot.margin = margin(5.5, 18, 5.5, 18)
  )

# save comparison figure
ggsave(
  filename = file.path(paths$figure_dir, "candidate_habitat_class_boxplots_2025.png"),
  plot = candidate_plot,
  device = ragg::agg_png,
  width = 12,
  height = 4,
  dpi = 300
)
