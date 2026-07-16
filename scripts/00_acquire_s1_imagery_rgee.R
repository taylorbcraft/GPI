# Build a Sentinel-1 temporal variability raster with rgee.
# The script filters the S1 GRD archive, smooths each scene in linear scale,
# computes consecutive log-ratio changes, and reduces those changes to a
# temporal standard deviation layer aligned to the Sentinel-2 predictor grid.

library(tidyverse)
library(sf)
library(terra)
library(rgee)
library(googledrive)

source("config.R")

# settings
s1_year <- Sys.getenv("S1_YEAR", unset = calibration_year)
start_date <- Sys.getenv("S1_START_DATE", unset = paste0(s1_year, "-03-31"))
end_date <- Sys.getenv("S1_END_DATE", unset = paste0(s1_year, "-07-17"))
orbit_pass <- Sys.getenv("S1_ORBIT_PASS", unset = "ASCENDING")
relative_orbit <- as.integer(Sys.getenv("S1_RELATIVE_ORBIT", unset = "88"))
polarization <- Sys.getenv("S1_POLARIZATION", unset = "VV")
instrument_mode <- Sys.getenv("S1_INSTRUMENT_MODE", unset = "IW")
boxcar_radius <- as.integer(Sys.getenv("S1_BOXCAR_RADIUS", unset = "3"))
export_scale <- as.integer(Sys.getenv("S1_EXPORT_SCALE", unset = "10"))
drive_folder <- "GEE"
study_area_source <- paths$field_geometry
output_dir <- file.path("scratch", "s1_rgee_imagery_test", "outputs", "data")
output_path <- file.path(
  paths$raster_dir,
  paste0(
    "S1_",
    polarization,
    "_LogRatio_StdDev_",
    str_to_title(str_to_lower(orbit_pass)),
    "_",
    s1_year,
    ".tif"
  )
)
metadata_path <- file.path(
  output_dir,
  paste0("s1_temporal_sd_metadata_", s1_year, ".csv")
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# rgee setup
ee_Initialize(drive = TRUE)

# study area
study_area <- st_read(study_area_source, quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(4326)

study_area_bbox <- study_area %>%
  st_union() %>%
  st_bbox()

study_area_ee <- ee$Geometry$Rectangle(
  coords = c(
    study_area_bbox[["xmin"]],
    study_area_bbox[["ymin"]],
    study_area_bbox[["xmax"]],
    study_area_bbox[["ymax"]]
  ),
  proj = "EPSG:4326",
  geodesic = FALSE
)

roi_extent_mask <- ee$Image$constant(1)$
  clip(study_area_ee)$
  selfMask()

# filter Sentinel-1 scenes
s1_collection <- ee$ImageCollection("COPERNICUS/S1_GRD")$
  filterBounds(study_area_ee)$
  filterDate(start_date, as.character(as.Date(end_date) + 1))$
  filter(ee$Filter$eq("orbitProperties_pass", orbit_pass))$
  filter(ee$Filter$inList("relativeOrbitNumber_start", list(relative_orbit)))$
  filter(ee$Filter$listContains("transmitterReceiverPolarisation", polarization))$
  filter(ee$Filter$eq("instrumentMode", instrument_mode))$
  select(polarization)$
  sort("system:time_start")

image_count <- s1_collection$size()$getInfo()

# spatial smoothing in linear power
boxcar <- ee$Kernel$square(
  radius = boxcar_radius,
  units = "pixels"
)

s1_linear_smoothed <- s1_collection$map(
  ee_utils_pyfunc(function(image) {
    linear_image <- ee$Image(10)$pow(image$divide(10))

    linear_image$
      convolve(boxcar)$
      rename("linear_smoothed")$
      copyProperties(image, image$propertyNames())
  })
)

# consecutive image pairs
image_list <- s1_linear_smoothed$toList(s1_linear_smoothed$size())
list_size <- s1_linear_smoothed$size()

log_ratio_list <- ee$List$sequence(0, list_size$subtract(2))$map(
  ee_utils_pyfunc(function(index) {
    i <- ee$Number(index)
    img1 <- ee$Image(image_list$get(i))
    img2 <- ee$Image(image_list$get(i$add(1)))

    img2$log10()$
      subtract(img1$log10())$
      rename("log_ratio")$
      set("system:time_start", img2$get("system:time_start"))
  })
)

log_ratio_collection <- ee$ImageCollection$fromImages(log_ratio_list)
log_ratio_count <- log_ratio_collection$size()$getInfo()

# reduce temporal variation
temporal_sd <- log_ratio_collection$
  reduce(ee$Reducer$stdDev())$
  updateMask(roi_extent_mask)$
  clip(study_area_ee)

# save run metadata
tibble(
  s1_year = s1_year,
  start_date = start_date,
  end_date = end_date,
  orbit_pass = orbit_pass,
  relative_orbit = relative_orbit,
  polarization = polarization,
  instrument_mode = instrument_mode,
  image_count = image_count,
  log_ratio_count = log_ratio_count,
  boxcar_radius = boxcar_radius,
  export_scale = export_scale
) %>%
  write_csv(metadata_path)

message("Images used: ", image_count)
message("Consecutive log ratios: ", log_ratio_count)
message("Orbit pass: ", orbit_pass)
message("Relative orbit: ", relative_orbit)

# export raster
download_path <- tempfile(
  pattern = paste0("s1_temporal_sd_", s1_year, "_raw_"),
  tmpdir = output_dir,
  fileext = ".tif"
)

ee_as_rast(
  image = temporal_sd,
  dsn = download_path,
  via = "drive",
  container = drive_folder,
  scale = export_scale,
  maxPixels = 1e13,
  timePrefix = FALSE,
  quiet = FALSE
)

# align to Sentinel-2 template grid
template_path <- path_calibration_raster("s2rep")
template_raster <- rast(template_path)
downloaded_raster <- rast(download_path)

if (!compareGeom(template_raster, downloaded_raster, stopOnError = FALSE)) {
  downloaded_raster <- project(downloaded_raster, template_raster, method = "bilinear")
}

writeRaster(downloaded_raster, output_path, overwrite = TRUE)
unlink(download_path)
