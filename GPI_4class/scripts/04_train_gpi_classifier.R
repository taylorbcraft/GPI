# Train the pixel-level random forest.
# Each pixel inherits its polygon's KNN-derived class label, and grouped
# cross-validation evaluates how well the model transfers to held-out polygons.

source("config.R")

library(tidyverse)
library(caret)
library(randomForest)
library(sf)
library(terra)
library(janitor)

rf_ntree <- 1000
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

# Test candidate mtry values
tune_mtry_repeated_cv <- function(data, predictors, mtry_values, seed = 123) {
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
    mtry = mtry_values,
    repeat_id = seq_len(cv_repeats),
    fold_id = seq_len(cv_folds)
  ) %>%
    mutate(
      metrics = pmap(
        list(mtry, repeat_id, fold_id),
        function(mtry, repeat_id, fold_id) {
          test_index <- resamples[[repeat_id]][[fold_id]]
          train_fold <- data[-test_index, ]
          test_fold <- data[test_index, ]

          rf_mod <- randomForest(
            x = train_fold %>% select(all_of(predictors)),
            y = train_fold$gpi_class,
            mtry = mtry,
            ntree = rf_ntree,
            importance = FALSE
          )

          pred <- predict(rf_mod, newdata = test_fold %>% select(all_of(predictors)))
          cm <- confusionMatrix(
            data = factor(pred, levels = gpi_class_levels),
            reference = factor(test_fold$gpi_class, levels = gpi_class_levels)
          )

          tibble(
            accuracy = unname(cm$overall["Accuracy"]),
            kappa = unname(cm$overall["Kappa"])
          )
        }
      )
    ) %>%
    unnest(metrics) %>%
    group_by(mtry) %>%
    summarise(
      Accuracy = mean(accuracy, na.rm = TRUE),
      Kappa = mean(kappa, na.rm = TRUE),
      AccuracySD = sd(accuracy, na.rm = TRUE),
      KappaSD = sd(kappa, na.rm = TRUE),
      n_resamples = dplyr::n(),
      .groups = "drop"
    ) %>%
    arrange(desc(Kappa), desc(Accuracy), mtry)
}

# Rerun grouped CV with selected mtry
predict_repeated_cv <- function(data, predictors, mtry, seed = 456) {
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
            mtry = mtry,
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

ensure_dirs(c(paths$validation_dir, paths$model_dir))

# Read KNN target
zone_labels <- read_csv(path_candidate_training(), show_col_types = FALSE)
target_var <- names(zone_labels) %>% str_subset("^gpi_class_") %>% dplyr::first()

zones <- st_read(paths$sampled_zone_geometry, quiet = TRUE) %>%
  clean_names() %>%
  mutate(
    polygon_id = as.character(polygon_id),
    meadow_id = str_remove(polygon_id, "_.*$")
  )

# Load calibration rasters
predictor_stack <- model_predictor_bands %>%
  set_names() %>%
  map(~ rast(path_calibration_raster(.x)))

if ("s2rep" %in% names(predictor_stack)) {
  predictor_stack$s2rep <- clamp(predictor_stack$s2rep, lower = 600, upper = 850, values = TRUE)
}

predictor_stack <- rast(predictor_stack)
zones <- st_transform(zones, crs(predictor_stack))

# Join polygon labels and extract pixel rows
labeled_zones <- zones %>%
  left_join(
    zone_labels %>%
      select(polygon_id, meadow_id, all_of(target_var)) %>%
      mutate(
        polygon_id = as.character(polygon_id),
        meadow_id = as.character(meadow_id),
        across(all_of(target_var), ~ factor(.x, levels = gpi_class_levels))
      ),
    by = c("polygon_id", "meadow_id")
  ) %>%
  filter(!is.na(.data[[target_var]]))

pixel_training_data <- terra::extract(
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
      select(zone_row, polygon_id, meadow_id, all_of(target_var)),
    by = "zone_row"
  ) %>%
  rename(gpi_class = all_of(target_var)) %>%
  drop_na(gpi_class, all_of(model_predictor_bands)) %>%
  mutate(
    gpi_class = factor(gpi_class, levels = gpi_class_levels),
    polygon_id = as.character(polygon_id),
    meadow_id = as.character(meadow_id)
  ) %>%
  select(cell_id, x, y, polygon_id, meadow_id, gpi_class, all_of(model_predictor_bands))

# Tune mtry
tuning_results <- tune_mtry_repeated_cv(
  data = pixel_training_data,
  predictors = model_predictor_bands,
  mtry_values = seq_along(model_predictor_bands),
  seed = 123
)

selected_mtry <- tuning_results %>%
  slice(1) %>%
  pull(mtry)

repeated_cv_predictions <- predict_repeated_cv(
  data = pixel_training_data,
  predictors = model_predictor_bands,
  mtry = selected_mtry,
  seed = 124
)

repeated_cv_eval <- summarise_confusion(
  pred = repeated_cv_predictions$prediction,
  ref = repeated_cv_predictions$reference,
  target_var = target_var
)

# Fit final forest
fitted_model <- randomForest(
  x = pixel_training_data %>% select(all_of(model_predictor_bands)),
  y = pixel_training_data$gpi_class,
  mtry = selected_mtry,
  ntree = rf_ntree,
  importance = TRUE
)

summary_tbl <- tibble(
  target = target_var,
  training_unit = "pixel",
  validation_group = "polygon_id",
  n_training_rows = nrow(pixel_training_data),
  n_training_polygons = n_distinct(pixel_training_data$polygon_id),
  selected_mtry = selected_mtry,
  ntree = rf_ntree,
  cv_folds = cv_folds,
  cv_repeats = cv_repeats,
  tuning_accuracy = tuning_results %>% slice(1) %>% pull(Accuracy),
  tuning_kappa = tuning_results %>% slice(1) %>% pull(Kappa),
  repeated_cv_predictions = nrow(repeated_cv_predictions),
  overall_accuracy = repeated_cv_eval$summary$overall_accuracy,
  kappa = repeated_cv_eval$summary$kappa
)

tuning_tbl <- tuning_results %>%
  as_tibble() %>%
  mutate(
    target = target_var,
    training_unit = "pixel",
    validation_group = "polygon_id",
    ntree = rf_ntree,
    cv_folds = cv_folds,
    cv_repeats = cv_repeats,
    selected = mtry == selected_mtry,
    .before = 1
  )

importance_tbl <- importance(fitted_model, type = 1) %>%
  as.data.frame() %>%
  rownames_to_column("predictor") %>%
  as_tibble() %>%
  rename(mean_decrease_accuracy = `MeanDecreaseAccuracy`) %>%
  mutate(target = target_var, .before = 1) %>%
  arrange(desc(mean_decrease_accuracy))

# Save diagnostics and model
write_csv(
  summary_tbl %>%
    mutate(
      best_target = target_var,
      compared_targets = target_var,
      predictors = paste(model_predictor_bands, collapse = ", "),
      .before = 1
    ),
  path_model_comparison()
)
write_csv(tuning_tbl, path_model_tuning())
write_csv(repeated_cv_eval$confusion, path_best_model_confusion())
write_csv(repeated_cv_eval$class_accuracy, path_best_model_class_accuracy())
write_csv(importance_tbl, path_best_model_variable_importance())
write_csv(
  summary_tbl %>%
    mutate(
      best_target = target_var,
      compared_targets = target_var,
      predictors = paste(model_predictor_bands, collapse = ", "),
      .before = 1
    ),
  path_best_model_metadata()
)

saveRDS(fitted_model, path_best_model())
