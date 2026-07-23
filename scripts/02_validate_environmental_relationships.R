# Validate how S2REP relates to the measured field variables.
# The script fits a joint GAM with a field-level random effect, formats the two
# manuscript GAM tables, and plots adjusted relationships with soil resistance,
# soil moisture, vegetation height, and plant richness.

source("config.R")

library(tidyverse)
library(mgcv)

ensure_dirs(c(paths$validation_dir, paths$figure_dir))
ensure_dirs("tables")

# field predictors
vars <- c(
  "soil_resistance",
  "soil_moisture",
  "vegetation_height",
  "plant_richness_mean"
)

# complete-case data
dat <- read_csv(path_anchor_training(), show_col_types = FALSE) %>%
  mutate(meadow_id = factor(meadow_id)) %>%
  select(meadow_id, s2rep, all_of(vars)) %>%
  drop_na()

# predictor labels
variable_labels <- c(
  soil_resistance = "soil resistance",
  soil_moisture = "soil moisture",
  vegetation_height = "vegetation height",
  plant_richness_mean = "plant richness"
)

# fit joint gam
gam_fit <- gam(
  s2rep ~
    s(soil_resistance, bs = "cs") +
    s(soil_moisture, bs = "cs") +
    s(vegetation_height, bs = "cs") +
    s(plant_richness_mean, bs = "cs") +
    s(meadow_id, bs = "re"),
  data = dat,
  method = "REML",
  select = TRUE
)

gam_fit_summary <- summary(gam_fit)
preds <- predict(gam_fit, type = "response")

# model fit summary
model_summary <- tibble(
  model = "gam",
  description = "Joint GAM with meadow random effect",
  n_obs = nobs(gam_fit),
  n_meadows = nlevels(dat$meadow_id),
  aic = AIC(gam_fit),
  bic = BIC(gam_fit),
  rmse = sqrt(mean((dat$s2rep - preds)^2)),
  r_squared = gam_fit_summary$r.sq,
  adj_r_squared = gam_fit_summary$r.sq,
  dev_explained = gam_fit_summary$dev.expl
)

# parametric and smooth terms
parametric_terms <- gam_fit_summary$p.table %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  as_tibble() %>%
  transmute(
    model = "gam",
    term_type = "parametric",
    term,
    estimate = Estimate,
    std_error = `Std. Error`,
    statistic = `t value`,
    p_value = `Pr(>|t|)`,
    edf = NA_real_,
    ref_df = NA_real_
  )

smooth_terms <- gam_fit_summary$s.table %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  as_tibble() %>%
  transmute(
    model = "gam",
    term_type = "smooth",
    term,
    estimate = NA_real_,
    std_error = NA_real_,
    statistic = F,
    p_value = `p-value`,
    edf = edf,
    ref_df = Ref.df
  )

model_coefficients <- bind_rows(parametric_terms, smooth_terms)

# full summary output
writeLines(
  c(
    paste0("Complete cases used: ", nrow(dat)),
    paste0("Unique meadows: ", nlevels(dat$meadow_id)),
    "",
    "Model summary",
    capture.output(print(model_summary)),
    "",
    "Joint generalized additive model with meadow random effect",
    capture.output(print(gam_fit_summary))
  ),
  file.path(paths$validation_dir, "environmental_validation_full_summary_2025.txt")
)

# manuscript table 1
model_summary %>%
  transmute(
    `n obs` = n_obs,
    `n fields` = n_meadows,
    AIC = round(aic, 3),
    BIC = round(bic, 3),
    RMSE = round(rmse, 3),
    `adj. R²` = round(adj_r_squared, 3),
    `dev. explained` = round(dev_explained, 3)
  ) %>%
  write_csv(file.path("tables", "table_1_gam_fit_statistics_2025.csv"))

# manuscript table 2
model_coefficients %>%
  transmute(
    `term type` = term_type,
    term = str_replace(term, "meadow_id", "field_id"),
    estimate = if_else(
      is.na(estimate),
      NA_character_,
      formatC(estimate, format = "f", digits = 3)
    ),
    statistic = case_when(
      statistic >= 1000 | statistic < 0.001 ~ formatC(statistic, format = "e", digits = 3),
      TRUE ~ formatC(statistic, format = "f", digits = 3)
    ),
    `p value` = case_when(
      p_value < 0.001 ~ "< 0.001",
      TRUE ~ formatC(p_value, format = "f", digits = 3)
    ),
    edf = if_else(
      is.na(edf),
      NA_character_,
      formatC(edf, format = "f", digits = 3)
    ),
    `ref. df` = if_else(
      is.na(ref_df),
      NA_character_,
      formatC(ref_df, format = "f", digits = 0)
    )
  ) %>%
  write_csv(
    file.path("tables", "table_2_gam_terms_2025.csv"),
    na = ""
  )

# effect plot data
partial_effect_data <- map_dfr(
  vars,
  function(var) {
    grid_vals <- seq(min(dat[[var]]), max(dat[[var]]), length.out = 200)

    newdata <- tibble(
      soil_resistance = median(dat$soil_resistance, na.rm = TRUE),
      soil_moisture = median(dat$soil_moisture, na.rm = TRUE),
      vegetation_height = median(dat$vegetation_height, na.rm = TRUE),
      plant_richness_mean = median(dat$plant_richness_mean, na.rm = TRUE),
      meadow_id = factor(levels(dat$meadow_id)[1], levels = levels(dat$meadow_id))
    ) %>%
      slice(rep(1, length(grid_vals))) %>%
      mutate(!!var := grid_vals)

    preds <- predict(
      gam_fit,
      newdata = newdata,
      type = "response",
      se.fit = TRUE,
      exclude = "s(meadow_id)"
    )

    tibble(
      predictor = var,
      x = grid_vals,
      fit = as.numeric(preds$fit),
      lower = fit - 1.96 * as.numeric(preds$se.fit),
      upper = fit + 1.96 * as.numeric(preds$se.fit)
    )
  }
) %>%
  mutate(
    predictor_label = recode(predictor, !!!variable_labels),
    predictor_label = factor(predictor_label, levels = unname(variable_labels[vars]))
  )

point_data <- dat %>%
  pivot_longer(cols = all_of(vars), names_to = "predictor", values_to = "x") %>%
  mutate(
    predictor_label = recode(predictor, !!!variable_labels),
    predictor_label = factor(predictor_label, levels = unname(variable_labels[vars]))
  )

# adjusted effect plot
combined_plot <- ggplot(partial_effect_data, aes(x = x, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey75", alpha = 0.7) +
  geom_line(color = "grey20", linewidth = 1) +
  geom_point(
    data = point_data,
    aes(x = x, y = s2rep),
    inherit.aes = FALSE,
    color = "grey35",
    alpha = 0.45,
    size = 1.5
  ) +
  facet_wrap(~ predictor_label, scales = "free_x") +
  labs(
    x = NULL,
    y = "Predicted S2REP"
  ) +
  theme_minimal(base_family = "Avenir Next", base_size = 16) +
  theme(
    text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.4),
    axis.text = element_text(color = "grey35"),
    axis.title = element_text(color = "grey20"),
    strip.text = element_text(face = "bold", size = 18, color = "grey20")
  )

# save validation plot
ggsave(
  filename = file.path(paths$figure_dir, "environmental_validation_plots_2025.png"),
  plot = combined_plot,
  device = ragg::agg_png,
  width = 10,
  height = 8,
  dpi = 300
)
