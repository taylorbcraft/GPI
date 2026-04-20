# Train and select the random forest classifier for the KNN-derived GPI target.
# The script tunes mtry and evaluates the selected setting with repeated
# cross-validation, then refits the saved model with all complete training rows.
#
# Input: data/processed/training/candidate_gpi_training_data_<calibration_year>.csv
# Outputs: model RDS, tuning tables, confusion matrix, class accuracy, variable
# importance, model comparison, and selected-model metadata.

source("config.R")

library(tidyverse)
library(caret)
library(randomForest)

rf_ntree <- 1000
cv_folds <- 5
cv_repeats <- 20

# Tune mtry with repeated folds so one lucky split does not drive model choice.
tune_mtry_repeated_cv <- function(data, target_var, predictors, mtry_values, seed = 123) {
  set.seed(seed)

  resamples <- map(
    seq_len(cv_repeats),
    ~ createFolds(
      data[[target_var]],
      k = cv_folds,
      returnTrain = FALSE
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
            y = train_fold[[target_var]],
            mtry = mtry,
            ntree = rf_ntree,
            importance = FALSE
          )

          pred <- predict(rf_mod, newdata = test_fold %>% select(all_of(predictors)))
          cm <- confusionMatrix(
            data = factor(pred, levels = gpi_class_levels),
            reference = factor(test_fold[[target_var]], levels = gpi_class_levels)
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

# Evaluate the selected setting with a fresh set of repeated folds.
predict_repeated_cv <- function(data, target_var, predictors, mtry, seed = 456) {
  set.seed(seed)

  resamples <- map(
    seq_len(cv_repeats),
    ~ createFolds(
      data[[target_var]],
      k = cv_folds,
      returnTrain = FALSE
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
            y = train_fold[[target_var]],
            mtry = mtry,
            ntree = rf_ntree,
            importance = FALSE
          )

          tibble(
            row_id = test_index,
            polygon_id = test_fold$polygon_id,
            meadow_id = test_fold$meadow_id,
            reference = test_fold[[target_var]],
            prediction = predict(rf_mod, newdata = test_fold %>% select(all_of(predictors)))
          )
        }
      )
    ) %>%
    unnest(predictions)
}

# Save overall and class-level metrics in a form that can be inspected later.
summarise_confusion <- function(pred, ref, target_var) {
  cm <- confusionMatrix(
    data = factor(pred, levels = gpi_class_levels),
    reference = factor(ref, levels = gpi_class_levels)
  )

  summary_tbl <- tibble(
    overall_accuracy = unname(cm$overall["Accuracy"]),
    kappa = unname(cm$overall["Kappa"])
  )

  class_acc_tbl <- tibble(
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
    )

  confusion_tbl <- as.data.frame(cm$table) %>%
    mutate(target = target_var, .before = 1)

  list(
    summary = summary_tbl,
    class_accuracy = class_acc_tbl,
    confusion = confusion_tbl
  )
}

ensure_dirs(c(paths$validation_dir, paths$model_dir))

require_file(path_candidate_training())

dat <- read_csv(path_candidate_training(), show_col_types = FALSE)

# Candidate target columns are created by script 03.
targets <- names(dat) %>%
  str_subset("^gpi_class_")

if (length(targets) == 0) {
  stop("No candidate gpi class columns found. Expected columns beginning with gpi_class_.", call. = FALSE)
}

# Only rows with complete model predictors can contribute to RF training.
dat_model <- dat %>%
  drop_na(all_of(model_predictor_bands))

# Fit one RF per candidate target; currently there is one configured target.
fit_rf <- function(data, target_var, predictors, seed = 123) {
  dat_sub <- data %>%
    select(polygon_id, meadow_id, all_of(target_var), all_of(predictors)) %>%
    drop_na(all_of(target_var)) %>%
    mutate(across(all_of(target_var), ~ factor(.x, levels = gpi_class_levels)))

  tuning_results <- tune_mtry_repeated_cv(
    data = dat_sub,
    target_var = target_var,
    predictors = predictors,
    mtry_values = seq_along(predictors),
    seed = seed
  )

  selected_mtry <- tuning_results %>%
    slice(1) %>%
    pull(mtry)

  repeated_cv_predictions <- predict_repeated_cv(
    data = dat_sub,
    target_var = target_var,
    predictors = predictors,
    mtry = selected_mtry,
    seed = seed + 1
  )

  repeated_cv_eval <- summarise_confusion(
    pred = repeated_cv_predictions$prediction,
    ref = repeated_cv_predictions$reference,
    target_var = target_var
  )

  summary_tbl <- tibble(
    target = target_var,
    n_training_rows = nrow(dat_sub),
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
      ntree = rf_ntree,
      cv_folds = cv_folds,
      cv_repeats = cv_repeats,
      selected = mtry == selected_mtry,
      .before = 1
    )

  final_model <- randomForest(
    x = dat_sub %>% select(all_of(predictors)),
    y = dat_sub[[target_var]],
    mtry = selected_mtry,
    ntree = rf_ntree,
    importance = TRUE
  )

  importance_tbl <- importance(final_model) %>%
    as.data.frame() %>%
    rownames_to_column("predictor") %>%
    mutate(target = target_var, .before = 1)

  list(
    model = final_model,
    summary = summary_tbl,
    tuning = tuning_tbl,
    class_accuracy = repeated_cv_eval$class_accuracy,
    importance = importance_tbl,
    confusion = repeated_cv_eval$confusion,
    training_data = dat_sub
  )
}

# Fit each candidate target and rank models by repeated-CV performance.
results <- map(
  targets,
  ~ fit_rf(
    data = dat_model,
    target_var = .x,
    predictors = model_predictor_bands,
    seed = 123
  )
)

names(results) <- targets

model_comparison <- map_dfr(results, "summary") %>%
  arrange(desc(kappa), desc(overall_accuracy))

model_tuning <- map_dfr(results, "tuning") %>%
  arrange(target, mtry)

best_target <- model_comparison %>%
  slice(1) %>%
  pull(target)

best_result <- results[[best_target]]
best_model <- best_result$model

# Metadata links the saved model to its target, predictors, and CV estimate.
best_model_metadata <- tibble(
  best_target = best_target,
  compared_targets = paste(targets, collapse = ", "),
  predictors = paste(model_predictor_bands, collapse = ", "),
  n_training_rows = nrow(best_result$training_data),
  selected_mtry = model_comparison %>% slice(1) %>% pull(selected_mtry),
  ntree = model_comparison %>% slice(1) %>% pull(ntree),
  cv_folds = cv_folds,
  cv_repeats = cv_repeats,
  tuning_accuracy = model_comparison %>% slice(1) %>% pull(tuning_accuracy),
  tuning_kappa = model_comparison %>% slice(1) %>% pull(tuning_kappa),
  repeated_cv_predictions = model_comparison %>% slice(1) %>% pull(repeated_cv_predictions),
  overall_accuracy = model_comparison %>% slice(1) %>% pull(overall_accuracy),
  kappa = model_comparison %>% slice(1) %>% pull(kappa)
)

write_csv(model_comparison, path_model_comparison())
write_csv(model_tuning, path_model_tuning())
write_csv(best_result$confusion, path_best_model_confusion())
write_csv(best_result$class_accuracy, path_best_model_class_accuracy())
write_csv(best_result$importance, path_best_model_variable_importance())
write_csv(best_model_metadata, path_best_model_metadata())
saveRDS(best_model, path_best_model())

print(model_comparison)
print(best_model_metadata)
