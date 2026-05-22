# Check whether the calibration table reflects sensible ecological structure.
# These plots and summaries are diagnostics only and do not feed into the classifier.

source("config.R")

library(tidyverse)
library(mgcv)
library(broom)
library(patchwork)

ensure_dirs(c(paths$validation_dir, paths$figure_dir))

# Read joined training table
dat <- read_csv(path_anchor_training(), show_col_types = FALSE)

vars <- c(
  "soil_resistance",
  "soil_moisture",
  "vegetation_height",
  "plant_richness_mean"
)

fit_models <- function(var) {
  f_lm <- as.formula(paste("s2rep ~", var))
  f_gam <- as.formula(paste("s2rep ~ s(", var, ")", sep = ""))

  mod_lm <- lm(f_lm, data = dat)
  mod_gam <- gam(f_gam, data = dat, method = "REML")

  tibble(
    predictor = var,
    aic_lm = AIC(mod_lm),
    aic_gam = AIC(mod_gam),
    delta_aic = abs(AIC(mod_lm) - AIC(mod_gam)),
    best_model = if ((AIC(mod_lm) - AIC(mod_gam)) >= 2) "gam" else "lm",
    r_squared_lm = summary(mod_lm)$r.squared,
    adj_r_squared_gam = summary(mod_gam)$r.sq
  )
}

model_summary <- map_dfr(vars, fit_models)
write_csv(model_summary, path_environmental_validation_summary())

# Use stronger model form in plots
plot_relationship <- function(var, xlab) {
  best_model <- model_summary %>%
    filter(predictor == var) %>%
    pull(best_model)

  p <- ggplot(dat, aes(x = .data[[var]], y = s2rep)) +
    geom_point(size = 2, alpha = 0.8) +
    labs(x = xlab, y = "s2rep") +
    theme_bw()

  if (best_model == "lm") {
    p + geom_smooth(method = "lm", se = TRUE)
  } else {
    p + geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = TRUE)
  }
}

combined_plot <- (
  plot_relationship("soil_resistance", "soil resistance") +
    plot_relationship("soil_moisture", "soil moisture")
) / (
  plot_relationship("vegetation_height", "vegetation height") +
    plot_relationship("plant_richness_mean", "plant richness")
)

ggsave(
  filename = path_environmental_validation_plot(),
  plot = combined_plot,
  width = 10,
  height = 8,
  dpi = 300
)

combined_plot
