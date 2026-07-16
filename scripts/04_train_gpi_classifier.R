# Train pixel-level random forests for each candidate target.
# Compare grouped cross-validation performance and save the best model outputs.

source("config.R")

library(tidyverse)
library(caret)
library(randomForest)
library(sf)
library(terra)
library(janitor)

rf_ntree <- 1000
rf_ntree_selection <- 250
rf_mtry <- 1
cv_folds <- 5
cv_repeats <- 20

if (!model_training_unit %in% c("pixel", "field")) {
  stop("model_training_unit must be 'pixel' or 'field'.")
}

if (length(candidate_model_predictor_bands) == 0) {
  stop("candidate_model_predictor_bands must contain at least one band.")
}

if (is.na(max_model_predictors) || max_model_predictors < 1) {
  stop("max_model_predictors must be at least 1.")
}

# polygon-grouped folds
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

# repeated stratified folds
make_stratified_folds <- function(labels, k, seed) {
  set.seed(seed)

  row_tbl <- tibble(
    row_id = seq_along(labels),
    label = as.character(labels)
  ) %>%
    arrange(label, row_id)

  fold_rows <- vector("list", k)

  for (label in unique(row_tbl$label)) {
    label_rows <- row_tbl %>%
      filter(label == !!label) %>%
      pull(row_id)

    label_rows <- sample(label_rows)
    assignments <- rep(seq_len(k), length.out = length(label_rows))

    for (fold_id in seq_len(k)) {
      fold_rows[[fold_id]] <- c(
        fold_rows[[fold_id]],
        label_rows[assignments == fold_id]
      )
    }
  }

  fold_rows
}

# repeated grouped cv
predict_repeated_cv_pixel <- function(data, predictors, seed = 456, cv_ntree = rf_ntree_selection) {
  model_mtry <- min(rf_mtry, length(predictors))

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
            mtry = model_mtry,
            ntree = cv_ntree,
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

predict_repeated_cv_field <- function(data, predictors, seed = 456, cv_ntree = rf_ntree_selection) {
  model_mtry <- min(rf_mtry, length(predictors))

  resamples <- map(
    seq_len(cv_repeats),
    ~ make_stratified_folds(
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
            mtry = model_mtry,
            ntree = cv_ntree,
            importance = FALSE
          )

          tibble(
            row_id = test_index,
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

# confusion-matrix outputs
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

# choose best candidate
choose_best_result <- function(results) {
  summary_tbl <- map_dfr(results, "summary") %>%
    mutate(result_row = row_number())

  best_row <- summary_tbl %>%
    arrange(
      desc(kappa),
      desc(overall_accuracy),
      n_predictors,
      predictors,
      target
    ) %>%
    slice(1) %>%
    pull(result_row)

  results[[best_row]]
}

# fit one target-predictor candidate
fit_target_model <- function(target_var, training_source, predictors, selection_stage,
                             cv_ntree = rf_ntree_selection, final_ntree = rf_ntree) {
  model_mtry <- min(rf_mtry, length(predictors))
  predictor_string <- paste(predictors, collapse = ", ")
  predictor_count <- length(predictors)

  if (model_training_unit == "pixel") {
    training_data <- training_source %>%
      transmute(
        cell_id,
        x,
        y,
        polygon_id,
        meadow_id,
        gpi_class = factor(.data[[target_var]], levels = gpi_class_levels),
        across(all_of(predictors))
      ) %>%
      drop_na(gpi_class, all_of(predictors))

    repeated_cv_predictions <- predict_repeated_cv_pixel(
      data = training_data,
      predictors = predictors,
      seed = 124,
      cv_ntree = cv_ntree
    )

    validation_group <- "polygon_id"
  } else {
    training_data <- training_source %>%
      transmute(
        polygon_id = as.character(polygon_id),
        meadow_id = as.character(meadow_id),
        gpi_class = factor(.data[[target_var]], levels = gpi_class_levels),
        across(all_of(predictors))
      ) %>%
      drop_na(gpi_class, all_of(predictors))

    repeated_cv_predictions <- predict_repeated_cv_field(
      data = training_data,
      predictors = predictors,
      seed = 124,
      cv_ntree = cv_ntree
    )

    validation_group <- "stratified_polygon_rows"
  }

  repeated_cv_eval <- summarise_confusion(
    pred = repeated_cv_predictions$prediction,
    ref = repeated_cv_predictions$reference,
    target_var = target_var
  )

  fitted_model <- randomForest(
    x = training_data %>% select(all_of(predictors)),
    y = training_data$gpi_class,
    mtry = model_mtry,
    ntree = final_ntree,
    importance = TRUE
  )

  summary_tbl <- tibble(
    target = target_var,
    selection_stage = selection_stage,
    predictors = predictor_string,
    n_predictors = predictor_count,
    training_unit = model_training_unit,
    validation_group = validation_group,
    n_training_rows = nrow(training_data),
    n_training_polygons = n_distinct(training_data$polygon_id),
    selected_mtry = model_mtry,
    ntree = final_ntree,
    ntree_selection = cv_ntree,
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
    mutate(
      target = target_var,
      predictors = predictor_string,
      .before = 1
    ) %>%
    arrange(desc(mean_decrease_accuracy))

  list(
    target = target_var,
    predictors = predictors,
    model = fitted_model,
    summary = summary_tbl,
    confusion = repeated_cv_eval$confusion,
    class_accuracy = repeated_cv_eval$class_accuracy,
    importance = importance_tbl
  )
}

# forward predictor search
run_forward_selection <- function(target_var, training_source) {
  results <- list()
  selected_predictors <- character()
  max_steps <- min(max_model_predictors, length(candidate_model_predictor_bands))

  for (selection_stage in seq_len(max_steps)) {
    if (selection_stage == 1) {
      predictor_sets <- map(candidate_model_predictor_bands, c)
    } else {
      remaining_predictors <- setdiff(candidate_model_predictor_bands, selected_predictors)

      if (length(remaining_predictors) == 0) {
        break
      }

      predictor_sets <- map(remaining_predictors, ~ c(selected_predictors, .x))
    }

    step_results <- map(
      predictor_sets,
      ~ fit_target_model(
        target_var = target_var,
        training_source = training_source,
        predictors = .x,
        selection_stage = selection_stage
      )
    )

    results <- c(results, step_results)
    selected_predictors <- choose_best_result(step_results)$predictors
  }

  results
}

ensure_dirs(c(paths$validation_dir, paths$model_dir))

# candidate target labels
zone_labels <- read_csv(path_candidate_training(), show_col_types = FALSE)
target_vars <- names(zone_labels) %>% str_subset("^gpi_class_")
zone_labels <- zone_labels %>%
  mutate(
    polygon_id = as.character(polygon_id),
    meadow_id = as.character(meadow_id),
    across(all_of(target_vars), ~ factor(.x, levels = gpi_class_levels))
  )

if (model_training_unit == "pixel") {
  zones <- st_read(paths$sampled_zone_geometry, quiet = TRUE) %>%
    clean_names() %>%
    mutate(
      polygon_id = as.character(polygon_id),
      meadow_id = str_remove(polygon_id, "_.*$")
    )

  predictor_stack <- load_predictor_stack(
    band_names = candidate_model_predictor_bands,
    image_date = calibration_image_date
  )
  zones <- st_transform(zones, crs(predictor_stack))

  labeled_zones <- zones %>%
    left_join(
      zone_labels %>%
        select(polygon_id, meadow_id, all_of(target_vars)),
      by = c("polygon_id", "meadow_id")
    )

  training_source <- terra::extract(
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
} else {
  anchor_training <- read_csv(path_anchor_training(), show_col_types = FALSE) %>%
    mutate(
      polygon_id = as.character(polygon_id),
      meadow_id = as.character(meadow_id)
    )

  training_source <- zone_labels %>%
    left_join(
      anchor_training %>%
        select(polygon_id, meadow_id, all_of(candidate_model_predictor_bands)),
      by = c("polygon_id", "meadow_id")
    )
}

# fit all candidate targets and predictor sets
target_results <- map(target_vars, ~ run_forward_selection(.x, training_source))
all_results <- flatten(target_results)
target_summary <- map_dfr(all_results, "summary")

best_result <- choose_best_result(all_results)
best_final_result <- fit_target_model(
  target_var = best_result$target,
  training_source = training_source,
  predictors = best_result$predictors,
  selection_stage = best_result$summary$selection_stage[[1]],
  cv_ntree = rf_ntree,
  final_ntree = rf_ntree
)

best_target <- best_final_result$target
compared_targets <- paste(target_vars, collapse = ", ")
compared_predictors <- paste(candidate_model_predictor_bands, collapse = ", ")
best_predictors <- paste(best_final_result$predictors, collapse = ", ")

# save diagnostics and model
write_csv(
  target_summary %>%
    mutate(
      best_target = best_target,
      best_predictors = best_predictors,
      compared_targets = compared_targets,
      compared_predictors = compared_predictors,
      .before = 1
    ),
  path_model_comparison()
)
write_csv(best_final_result$confusion, path_best_model_confusion())
write_csv(best_final_result$class_accuracy, path_best_model_class_accuracy())
write_csv(best_final_result$importance, path_best_model_variable_importance())
write_csv(
  best_final_result$summary %>%
    mutate(
      best_target = best_target,
      best_predictors = best_predictors,
      compared_targets = compared_targets,
      compared_predictors = compared_predictors,
      .before = 1
    ),
  path_best_model_metadata()
)

saveRDS(best_final_result$model, path_best_model())
