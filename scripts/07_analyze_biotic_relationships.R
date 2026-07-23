# Join the 2025 field classes to the godwit and insect datasets.
# The script summarizes godwit occupancy, occupied-field density, and habitat
# selection, links insect biomass and carabid richness to field classes, and
# creates the biodiversity figures.

source("config.R")

library(tidyverse)
library(sf)
library(ggplot2)
library(patchwork)

plot_by_class <- function(data, value_col, y_label) {
  ggplot(data, aes(x = habitat_class, y = .data[[value_col]], fill = habitat_class)) +
    geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85, linewidth = 0.8) +
    geom_jitter(width = 0.12, alpha = 0.45, size = 1.4, color = "grey20") +
    scale_fill_manual(
      values = habitat_class_palette,
      breaks = habitat_class_levels,
      drop = FALSE
    ) +
    labs(
      x = "Habitat quality",
      y = y_label
    ) +
    theme_minimal(base_family = "Avenir Next", base_size = 12) +
    theme(
      text = element_text(face = "bold"),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.4),
      axis.text = element_text(color = "grey35"),
      axis.title = element_text(color = "grey20")
    )
}

ensure_dirs(paths$figure_dir)

# field classes
field_classes <- st_read(path_field_map(prediction_year), quiet = TRUE) %>%
  st_drop_geometry() %>%
  transmute(
    meadow_id = as.character(field_id),
    s2rep_mean,
    habitat_class = factor(habitat_class, levels = habitat_class_levels),
    class_rank = match(habitat_class, habitat_class_levels)
  )

# godwit study fields
polders <- as.character(c(
  14, 16, 17, 66, 149, 154, 155, 156, 5, 7, 24, 116, 114, 583, 6, 12, 25,
  100, 147, 580, 79, 172, 258, 563, 564, 565, 566, 567, 568, 569, 570, 571,
  572, 573, 574, 575, 576, 577, 578, 115, 142, 143, 144, 145, 146, 158, 579,
  581, 582, 129, 162, 69, 98, 127, 128, 130, 157, 159, 160, 161, 164, 70,
  131, 133, 134, 137, 138, 150, 163
))

polder_pattern <- paste0("^(", paste(polders, collapse = "|"), ")\\d{3}$")

# godwit survey counts
godwit_counts <- read_csv(
  file.path("data", "raw", "biotic", "pc_btg.csv"),
  show_col_types = FALSE
) %>%
  filter(Species == "BTG_terr") %>%
  transmute(
    meadow_id = as.character(MeadowID),
    round = factor(as.character(Ronde), levels = c("1", "2", "3")),
    count = as.numeric(Number_Counted)
  )

# complete field-round grid
godwit_rounds <- read_csv(
  file.path("data", "raw", "biotic", "meadow_ha.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    meadow_id = as.character(meadowID),
    area_ha
  ) %>%
  filter(str_detect(meadow_id, polder_pattern)) %>%
  crossing(round = factor(c("1", "2", "3"), levels = c("1", "2", "3"))) %>%
  left_join(field_classes, by = "meadow_id") %>%
  filter(!is.na(habitat_class)) %>%
  left_join(godwit_counts, by = c("meadow_id", "round")) %>%
  mutate(count = replace_na(count, 0))

# field-level godwit responses
godwit_fields <- godwit_rounds %>%
  group_by(meadow_id) %>%
  summarise(
    area_ha = first(area_ha),
    habitat_class = first(habitat_class),
    class_rank = first(class_rank),
    pc_btg = mean(count, na.rm = TRUE),
    godwit_density_ha = mean(count / area_ha, na.rm = TRUE),
    occupied = any(count > 0),
    .groups = "drop"
  )

# habitat selection ratios
godwit_selection <- godwit_fields %>%
  group_by(habitat_class, class_rank) %>%
  summarise(
    total_used = sum(pc_btg, na.rm = TRUE),
    total_avail = sum(area_ha, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    prop_used = total_used / sum(total_used),
    prop_avail = total_avail / sum(total_avail),
    wi = prop_used / prop_avail
  ) %>%
  arrange(class_rank) %>%
  select(-class_rank)

# occupancy summary
godwit_occupancy <- godwit_fields %>%
  group_by(habitat_class, class_rank) %>%
  summarise(
    n_fields = n(),
    occupied_fields = sum(occupied),
    occupancy_rate = occupied_fields / n_fields,
    conf_low = prop.test(occupied_fields, n_fields)$conf.int[[1]],
    conf_high = prop.test(occupied_fields, n_fields)$conf.int[[2]],
    .groups = "drop"
  ) %>%
  arrange(class_rank) %>%
  select(-class_rank)

# occupied-field density summary
godwit_occupied_density <- godwit_fields %>%
  filter(occupied) %>%
  group_by(habitat_class) %>%
  summarise(
    n_occupied = n(),
    mean_density_ha = mean(godwit_density_ha),
    se_density_ha = sd(godwit_density_ha) / sqrt(n_occupied),
    conf_low = pmax(mean_density_ha - 1.96 * se_density_ha, 0),
    conf_high = mean_density_ha + 1.96 * se_density_ha,
    .groups = "drop"
  )

# insect biomass
insect_presence <- read_csv(
  file.path("data", "raw", "biotic", "insect_data.csv"),
  show_col_types = FALSE
) %>%
  distinct(MeadowID) %>%
  transmute(meadow_id = as.character(MeadowID))

insect_biomass <- read.csv(file.path("data", "raw", "biotic", "insect_data_mass.csv")) %>%
  as_tibble() %>%
  transmute(
    meadow_id = as.character(MeadowID),
    insect_biomass_g = Gewicht.monster..g.
  ) %>%
  filter(!is.na(insect_biomass_g)) %>%
  group_by(meadow_id) %>%
  summarise(
    insect_biomass_g = mean(insect_biomass_g, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  right_join(insect_presence, by = "meadow_id") %>%
  mutate(
    insect_biomass_g = replace_na(insect_biomass_g, 0),
    n_samples = replace_na(n_samples, 0L)
  ) %>%
  left_join(field_classes, by = "meadow_id") %>%
  filter(!is.na(habitat_class))

# carabid richness
carabid_richness <- read.csv(file.path("data", "raw", "biotic", "insect_data_coleoptera.csv")) %>%
  as_tibble() %>%
  transmute(
    meadow_id = as.character(MeadowID),
    carabid_group = Familie.Geslacht
  ) %>%
  filter(carabid_group %in% c(
    "Acupalpus", "Amara", "Anisodactilus", "Bembidion", "Carabus", "Clivina",
    "Dyschirius", "Harpalus", "Loricera", "Oodes", "Poecilus", "Pterostichus",
    "Stomis"
  )) %>%
  distinct(meadow_id, carabid_group) %>%
  count(meadow_id, name = "carabid_family_richness") %>%
  right_join(insect_presence, by = "meadow_id") %>%
  mutate(carabid_family_richness = replace_na(carabid_family_richness, 0L)) %>%
  left_join(field_classes, by = "meadow_id") %>%
  filter(!is.na(habitat_class))

# occupancy and density panels
godwit_plot <- ((
  ggplot(godwit_occupancy, aes(x = habitat_class, y = occupancy_rate, fill = habitat_class)) +
    geom_col(width = 0.65, alpha = 0.9) +
    geom_errorbar(aes(ymin = conf_low, ymax = conf_high), width = 0.15) +
    scale_fill_manual(
      values = habitat_class_palette,
      breaks = habitat_class_levels,
      drop = FALSE
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      x = "Habitat quality",
      y = "Occupancy rate"
    ) +
    theme_minimal(base_family = "Avenir Next", base_size = 12)
) | (
  ggplot(godwit_occupied_density, aes(x = habitat_class, y = mean_density_ha, color = habitat_class)) +
    geom_pointrange(
      aes(ymin = conf_low, ymax = conf_high),
      linewidth = 0.7,
      size = 0.5
    ) +
    scale_color_manual(
      values = habitat_class_palette,
      breaks = habitat_class_levels,
      drop = FALSE
    ) +
    labs(
      x = "Habitat quality",
      y = "Mean territorial godwits ha⁻¹"
    ) +
    theme_minimal(base_family = "Avenir Next", base_size = 12)
)) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(
    text = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.4),
    axis.text = element_text(color = "grey35"),
    axis.title = element_text(color = "grey20"),
    plot.tag = element_text(
      family = "Avenir Next",
      face = "bold",
      size = 11,
      color = "grey20"
    )
  )

# habitat selection figure
godwit_selection_plot <- ggplot(
  godwit_selection,
  aes(x = habitat_class, y = wi, fill = habitat_class)
) +
  geom_col(width = 0.9, color = "black", linewidth = 0.4) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.5) +
  scale_fill_manual(
    values = habitat_class_palette,
    breaks = habitat_class_levels,
    drop = FALSE
  ) +
  labs(
    x = "Habitat quality",
    y = "Manly habitat selection ratio"
  ) +
  theme_minimal(base_family = "Avenir Next", base_size = 12) +
  theme(
    text = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.4),
    axis.text = element_text(color = "grey35"),
    axis.title = element_text(color = "grey20")
  )

# insect response panels
insect_plot <- ((
  plot_by_class(insect_biomass, "insect_biomass_g", "Aerial insect biomass (g)") +
    ggtitle("Aerial insect biomass")
) | (
  plot_by_class(carabid_richness, "carabid_family_richness", "Carabid family richness") +
  ggtitle("Carabid family richness")
)) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(
    plot.title = element_text(
      family = "Avenir Next",
      face = "bold",
      color = "grey20"
    ),
    plot.tag = element_text(
      family = "Avenir Next",
      face = "bold",
      size = 11,
      color = "grey20"
    )
  )

# save figures
ggsave(
  filename = file.path(paths$figure_dir, paste0("godwit_relationships_", prediction_year, ".png")),
  plot = godwit_plot,
  device = ragg::agg_png,
  width = 6,
  height = 4,
  dpi = 300
)

ggsave(
  filename = file.path(paths$figure_dir, paste0("godwit_selection_ratio_", prediction_year, ".png")),
  plot = godwit_selection_plot,
  device = ragg::agg_png,
  width = 5.5,
  height = 4,
  dpi = 300
)

ggsave(
  filename = file.path(paths$figure_dir, paste0("insect_relationships_", prediction_year, ".png")),
  plot = insect_plot,
  device = ragg::agg_png,
  width = 8,
  height = 4,
  dpi = 300
)
