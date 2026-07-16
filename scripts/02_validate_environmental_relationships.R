# Validate how S2REP relates to the measured field variables.
# Fit one joint GAM, write summary outputs, and save adjusted-effect plots.

source("config.R")

library(tidyverse)
library(mgcv)

ensure_dirs(c(paths$validation_dir, paths$figure_dir))

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
full_model_summaries <- c(
  paste0("Complete cases used: ", nrow(dat)),
  paste0("Unique meadows: ", nlevels(dat$meadow_id)),
  "",
  "Model summary",
  capture.output(print(model_summary)),
  "",
  "Joint generalized additive model with meadow random effect",
  capture.output(print(gam_fit_summary))
)

write_csv(model_summary, path_environmental_validation_summary())
write_csv(model_coefficients, path_environmental_validation_coefficients())
writeLines(full_model_summaries, path_environmental_validation_full_summary())

# adjusted effect curves
make_partial_effect_data <- function(var) {
  grid_vals <- seq(min(dat[[var]]), max(dat[[var]]), length.out = 200)
  base_vals <- set_names(vars) %>%
    map_dbl(~ median(dat[[.x]], na.rm = TRUE)) %>%
    as.list()

  newdata <- as_tibble(base_vals) %>%
    slice(rep(1, length(grid_vals))) %>%
    mutate(
      meadow_id = factor(levels(dat$meadow_id)[1], levels = levels(dat$meadow_id)),
      !!var := grid_vals
    )

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

# effect plot data
partial_effect_data <- map_dfr(vars, make_partial_effect_data) %>%
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
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#c7dbe6", alpha = 0.7) +
  geom_line(color = "#184e77", linewidth = 1) +
  geom_point(
    data = point_data,
    aes(x = x, y = s2rep),
    inherit.aes = FALSE,
    color = "black",
    alpha = 0.45,
    size = 1.5
  ) +
  facet_wrap(~ predictor_label, scales = "free_x") +
  labs(
    title = "S2REP relationships with soil and vegetation properties",
    subtitle = "Points are raw observed data; lines show fitted relationships holding other field variables at their median values",
    x = NULL,
    y = "Predicted S2REP"
  ) +
  theme_bw()

# save validation plot
ggsave(
  filename = path_environmental_validation_plot(),
  plot = combined_plot,
  width = 10,
  height = 8,
  dpi = 300
)
  
combined_plot
