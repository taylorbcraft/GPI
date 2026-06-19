# Train pixel-level random forests for each candidate target
# Each pixel inherits its polygon's candidate class label, and grouped
# cross-validation compares which target transfers best to held-out polygons

source("config.R")

library(tidyverse)
library(caret)
library(randomForest)
library(sf)
library(terra)
library(janitor)

rf_ntree <- 1000
rf_mtry <- 1
cv_folds <- 5
cv_repeats <- 20

# Keep whole polygons together across folds
make_grouped_folds <- function(groups, labels, k, seed) {
  set.seed(seed)

  group_tbl <- tibble(
    group = as.character(groups),
    label = as.character(labels)
  ) %>%
    distinct(group, label) %>%
    arrange(label, group)

  fold_groups <- vector("list", k)

  for (label in unique(group_tbl$label)) {
    label_groups <- group_tbl %>%
      filter(label == !!label) %>%
      pull(group)

    label_groups <- sample(label_groups)
    assignments <- rep(seq_len(k), length.out = length(label_groups))

    for (fold_id in seq_len(k)) {
      fold_groups[[fold_id]] <- c(
        fold_groups[[fold_id]],
        label_groups[assignments == fold_id]
      )
    }
  }

  map(fold_groups, ~ which(as.character(groups) %in% .x))
}

# Rerun grouped CV
predict_repeated_cv <- function(data, predictors, seed = 456) {
  resamples <- map(
    seq_len(cv_repeats),
    ~ make_grouped_folds(
      groups = data$polygon_id,
      labels = data$gpi_class,
      k = cv_folds,
      seed = seed + .x - 1
    )
  )

  crossing(
    repeat_id = seq_len(cv_repeats),
    fold_id = seq_len(cv_folds)
  ) %>%
    mutate(
      predictions = map2(
        repeat_id,
        fold_id,
        function(repeat_id, fold_id) {
          test_index <- resamples[[repeat_id]][[fold_id]]
          train_fold <- data[-test_index, ]
          test_fold <- data[test_index, ]

          rf_mod <- randomForest(
            x = train_fold %>% select(all_of(predictors)),
            y = train_fold$gpi_class,
            mtry = rf_mtry,
            ntree = rf_ntree,
            importance = FALSE
          )

          tibble(
            row_id = test_index,
            cell_id = test_fold$cell_id,
            polygon_id = test_fold$polygon_id,
            meadow_id = test_fold$meadow_id,
            reference = test_fold$gpi_class,
            prediction = predict(rf_mod, newdata = test_fold %>% select(all_of(predictors)))
          )
        }
      )
    ) %>%
    unnest(predictions)
}

# Turn confusion matrix into output tables
summarise_confusion <- function(pred, ref, target_var) {
  cm <- confusionMatrix(
    data = factor(pred, levels = gpi_class_levels),
    reference = factor(ref, levels = gpi_class_levels)
  )

  list(
    summary = tibble(
      overall_accuracy = unname(cm$overall["Accuracy"]),
      kappa = unname(cm$overall["Kappa"])
    ),
    class_accuracy = tibble(
      class = rownames(cm$byClass),
      sensitivity = cm$byClass[, "Sensitivity"],
      specificity = cm$byClass[, "Specificity"],
      pos_pred_value = cm$byClass[, "Pos Pred Value"],
      neg_pred_value = cm$byClass[, "Neg Pred Value"],
      balanced_accuracy = cm$byClass[, "Balanced Accuracy"]
    ) %>%
      rename(
        producers_accuracy = sensitivity,
        users_accuracy = pos_pred_value
      ) %>%
      mutate(
        target = target_var,
        class = str_remove(class, "^Class: "),
        .before = 1
      ),
    confusion = as.data.frame(cm$table) %>%
      mutate(target = target_var, .before = 1)
  )
}

# Fit and score one candidate target
fit_target_model <- function(target_var, pixel_labels) {
  pixel_training_data <- pixel_labels %>%
    transmute(
      cell_id,
      x,
      y,
      polygon_id,
      meadow_id,
      gpi_class = factor(.data[[target_var]], levels = gpi_class_levels),
      across(all_of(model_predictor_bands))
    ) %>%
    drop_na(gpi_class, all_of(model_predictor_bands))

  repeated_cv_predictions <- predict_repeated_cv(
    data = pixel_training_data,
    predictors = model_predictor_bands,
    seed = 124
  )

  repeated_cv_eval <- summarise_confusion(
    pred = repeated_cv_predictions$prediction,
    ref = repeated_cv_predictions$reference,
    target_var = target_var
  )

  fitted_model <- randomForest(
    x = pixel_training_data %>% select(all_of(model_predictor_bands)),
    y = pixel_training_data$gpi_class,
    mtry = rf_mtry,
    ntree = rf_ntree,
    importance = TRUE
  )

  summary_tbl <- tibble(
    target = target_var,
    training_unit = "pixel",
    validation_group = "polygon_id",
    n_training_rows = nrow(pixel_training_data),
    n_training_polygons = n_distinct(pixel_training_data$polygon_id),
    selected_mtry = rf_mtry,
    ntree = rf_ntree,
    cv_folds = cv_folds,
    cv_repeats = cv_repeats,
    repeated_cv_predictions = nrow(repeated_cv_predictions),
    overall_accuracy = repeated_cv_eval$summary$overall_accuracy,
    kappa = repeated_cv_eval$summary$kappa
  )

  importance_tbl <- importance(fitted_model, type = 1) %>%
    as.data.frame() %>%
    rownames_to_column("predictor") %>%
    as_tibble() %>%
    rename(mean_decrease_accuracy = `MeanDecreaseAccuracy`) %>%
    mutate(target = target_var, .before = 1) %>%
    arrange(desc(mean_decrease_accuracy))

  list(
    target = target_var,
    model = fitted_model,
    summary = summary_tbl,
    confusion = repeated_cv_eval$confusion,
    class_accuracy = repeated_cv_eval$class_accuracy,
    importance = importance_tbl
  )
}

ensure_dirs(c(paths$validation_dir, paths$model_dir))

# Read candidate targets
zone_labels <- read_csv(path_candidate_training(), show_col_types = FALSE)
target_vars <- names(zone_labels) %>% str_subset("^gpi_class_")

zones <- st_read(paths$sampled_zone_geometry, quiet = TRUE) %>%
  clean_names() %>%
  mutate(
    polygon_id = as.character(polygon_id),
    meadow_id = str_remove(polygon_id, "_.*$")
  )

# Load calibration rasters
predictor_stack <- load_predictor_stack(
  band_names = model_predictor_bands,
  image_date = calibration_image_date
)
zones <- st_transform(zones, crs(predictor_stack))

# Join polygon labels and extract pixel rows once
labeled_zones <- zones %>%
  left_join(
    zone_labels %>%
      select(polygon_id, meadow_id, all_of(target_vars)) %>%
      mutate(
        polygon_id = as.character(polygon_id),
        meadow_id = as.character(meadow_id),
        across(all_of(target_vars), ~ factor(.x, levels = gpi_class_levels))
      ),
    by = c("polygon_id", "meadow_id")
  )

pixel_labels <- terra::extract(
  x = predictor_stack,
  y = vect(labeled_zones),
  cells = TRUE,
  xy = TRUE
) %>%
  as_tibble() %>%
  rename(zone_row = ID, cell_id = cell) %>%
  left_join(
    st_drop_geometry(labeled_zones) %>%
      mutate(zone_row = row_number()) %>%
      select(zone_row, polygon_id, meadow_id, all_of(target_vars)),
    by = "zone_row"
  ) %>%
  mutate(
    polygon_id = as.character(polygon_id),
    meadow_id = as.character(meadow_id)
  )

# Fit all candidate targets
target_results <- map(target_vars, ~ fit_target_model(.x, pixel_labels))
target_summary <- map_dfr(target_results, "summary")

best_result <- target_results %>%
  set_names(target_vars) %>%
  pluck(
    target_summary %>%
      arrange(desc(kappa), desc(overall_accuracy), target) %>%
      slice(1) %>%
      pull(target)
  )

best_target <- best_result$target
compared_targets <- paste(target_vars, collapse = ", ")

# Save diagnostics and model
write_csv(
  target_summary %>%
    mutate(
      best_target = best_target,
      compared_targets = compared_targets,
      predictors = paste(model_predictor_bands, collapse = ", "),
      .before = 1
    ),
  path_model_comparison()
)
write_csv(best_result$confusion, path_best_model_confusion())
write_csv(best_result$class_accuracy, path_best_model_class_accuracy())
write_csv(best_result$importance, path_best_model_variable_importance())
write_csv(
  best_result$summary %>%
    mutate(
      best_target = best_target,
      compared_targets = compared_targets,
      predictors = paste(model_predictor_bands, collapse = ", "),
      .before = 1
    ),
  path_best_model_metadata()
)

saveRDS(best_result$model, path_best_model())
