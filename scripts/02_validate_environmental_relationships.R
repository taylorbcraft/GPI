# Validate whether the anchor table contains interpretable ecological signals.
# For each measured field variable, compare a linear model and a GAM for s2rep.
# The summary and plot are diagnostics; they do not feed directly into training.
#
# Input: data/processed/training/anchor_zone_training_data_<calibration_year>.csv
# Outputs: environmental validation CSV and figure.

source("config.R")

library(tidyverse)
library(mgcv)
library(broom)
library(patchwork)

ensure_dirs(c(paths$validation_dir, paths$figure_dir))

require_file(path_anchor_training())

dat <- read_csv(path_anchor_training(), show_col_types = FALSE)

# These diagnostics check whether S2REP tracks the measured field gradients.
vars <- c(
  "soil_resistance",
  "soil_moisture",
  "vegetation_height",
  "plant_richness_mean"
)

# Compare a simple linear response with a smooth response for each variable.
fit_models <- function(var) {
  f_lm <- as.formula(paste("s2rep ~", var))
  f_gam <- as.formula(paste("s2rep ~ s(", var, ")", sep = ""))

  mod_lm <- lm(f_lm, data = dat)
  mod_gam <- gam(f_gam, data = dat, method = "REML")

  aic_lm <- AIC(mod_lm)
  aic_gam <- AIC(mod_gam)
  best_model <- if ((aic_lm - aic_gam) >= 2) "gam" else "lm"

  tibble(
    predictor = var,
    aic_lm = aic_lm,
    aic_gam = aic_gam,
    delta_aic = abs(aic_lm - aic_gam),
    best_model = best_model,
    r_squared_lm = summary(mod_lm)$r.squared,
    adj_r_squared_gam = summary(mod_gam)$r.sq
  )
}

model_summary <- map_dfr(vars, fit_models)

write_csv(model_summary, path_environmental_validation_summary())

# Plot the relationship using whichever model form had stronger AIC support.
plot_relationship <- function(var, xlab) {
  row <- model_summary %>%
    filter(predictor == var)

  p <- ggplot(dat, aes(x = .data[[var]], y = s2rep)) +
    geom_point(size = 2, alpha = 0.8) +
    labs(x = xlab, y = "s2rep") +
    theme_bw()

  if (row$best_model == "lm") {
    p <- p + geom_smooth(method = "lm", se = TRUE)
  } else {
    p <- p + geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = TRUE)
  }

  p
}

p1 <- plot_relationship("soil_resistance", "soil resistance")
p2 <- plot_relationship("soil_moisture", "soil moisture")
p3 <- plot_relationship("vegetation_height", "vegetation height")
p4 <- plot_relationship("plant_richness_mean", "plant richness")
combined_plot <- (p1 + p2) / (p3 + p4)

ggsave(
  filename = path_environmental_validation_plot(),
  plot = combined_plot,
  width = 10,
  height = 8,
  dpi = 300
)

combined_plot
