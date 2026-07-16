library(tidyverse)
library(sf)
library(terra)
library(rgee)
library(lubridate)
library(googledrive)

# Acquire the project Sentinel-2 predictor rasters with rgee.
source("config.R")

# settings
composite_mode <- Sys.getenv("S2_COMPOSITE_MODE", unset = "monthly_median")
target_date <- Sys.getenv("S2_TARGET_DATE", unset = calibration_image_date)
fill_date <- Sys.getenv("S2_FILL_DATE", unset = "2025-04-06")
apply_smoothing <- Sys.getenv("S2_APPLY_SMOOTHING", unset = "true")
smoothing_radius <- as.integer(Sys.getenv("S2_SMOOTHING_RADIUS", unset = "3"))
mask_shadows <- Sys.getenv("S2_MASK_SHADOWS", unset = "false")

if (!composite_mode %in% c("monthly_median", "single_day")) {
  stop("S2_COMPOSITE_MODE must be 'monthly_median' or 'single_day'.")
}

if (!apply_smoothing %in% c("true", "false")) {
  stop("S2_APPLY_SMOOTHING must be 'true' or 'false'.")
}

if (!mask_shadows %in% c("true", "false")) {
  stop("S2_MASK_SHADOWS must be 'true' or 'false'.")
}

if (is.na(smoothing_radius) || smoothing_radius < 0) {
  stop("S2_SMOOTHING_RADIUS must be a non-negative integer.")
}

target_month <- as.Date(if (grepl("-median$", target_date)) {
  sub("-median$", "-01", target_date)
} else {
  target_date
})

start_date <- floor_date(target_month, unit = "month")
end_date <- ceiling_date(target_month, unit = "month") - days(1)
export_label <- if (composite_mode == "monthly_median") {
  paste0(format(start_date, "%Y-%m"), "-median")
} else {
  target_date
}
cloud_threshold <- 65
export_scale <- 20
drive_folder <- "GEE"
study_area_source <- paths$field_geometry
output_dir <- file.path("scratch", "s2_rgee_imagery_test", "outputs", "data")
export_bands <- c("s2rep", "msi")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# rgee setup
ee_Initialize(drive = TRUE)

# drive auth
if (interactive() && !drive_has_token()) {
  drive_auth()
}

if (!drive_has_token()) {
  stop(
    paste(
      "Google Drive auth is not available in this session.",
      "Run this script from an interactive R session after authenticating Drive,",
      "or call googledrive::drive_auth() first."
    )
  )
}

# study area
study_area <- st_read(study_area_source, quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(4326)

study_area_bbox <- study_area %>%
  st_union() %>%
  st_bbox() %>%
  unclass()

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

# source collections
s2 <- ee$ImageCollection("COPERNICUS/S2_SR_HARMONIZED")$
  filterBounds(study_area_ee)$
  filterDate(as.character(start_date), as.character(end_date + days(1)))

s2_clouds <- ee$ImageCollection("COPERNICUS/S2_CLOUD_PROBABILITY")$
  filterBounds(study_area_ee)$
  filterDate(as.character(start_date), as.character(end_date + days(1)))

joined <- ee$Join$saveFirst("cloud_mask")$apply(
  primary = s2,
  secondary = s2_clouds,
  condition = ee$Filter$equals(
    leftField = "system:index",
    rightField = "system:index"
  )
)

joined_collection <- ee$ImageCollection(joined)

# cloud mask
mask_clouds <- ee_utils_pyfunc(function(img) {
  cloud_prob <- ee$Image(img$get("cloud_mask"))$select("probability")
  cloud_mask <- cloud_prob$lt(cloud_threshold)
  combined_mask <- cloud_mask

  if (mask_shadows == "true") {
    shadow_mask <- img$select("SCL")$neq(3L)
    combined_mask <- cloud_mask$And(shadow_mask)
  }

  img$
    updateMask(combined_mask)$
    divide(10000)$
    copyProperties(img, img$propertyNames())
})

# predictor bands
add_indices <- ee_utils_pyfunc(function(img) {
  s2rep <- img$expression(
    "705 + 35 * ((((R + RE3) / 2) - RE1) / (RE2 - RE1))",
    list(
      R = img$select("B4"),
      RE1 = img$select("B5"),
      RE2 = img$select("B6"),
      RE3 = img$select("B7")
    )
  )$rename("s2rep")

  ndvi <- img$normalizedDifference(c("B8", "B4"))$rename("ndvi")
  ndwi <- img$normalizedDifference(c("B3", "B8"))$rename("ndwi")

  savi <- img$expression(
    "((NIR - RED) / (NIR + RED + L)) * (1 + L)",
    list(
      NIR = img$select("B8"),
      RED = img$select("B4"),
      L = 0.5
    )
  )$rename("savi")

  evi <- img$expression(
    "2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))",
    list(
      NIR = img$select("B8"),
      RED = img$select("B4"),
      BLUE = img$select("B2")
    )
  )$rename("evi")

  msi <- img$expression(
    "SWIR1 / NIR",
    list(
      SWIR1 = img$select("B11"),
      NIR = img$select("B8")
    )
  )$rename("msi")

  ndmi <- img$normalizedDifference(c("B8", "B11"))$rename("ndmi")
  mndwi <- img$normalizedDifference(c("B3", "B11"))$rename("mndwi")

  img$addBands(list(s2rep, ndvi, ndwi, savi, evi, msi, ndmi, mndwi))
})

processed <- joined_collection$
  map(mask_clouds)$
  map(add_indices)

make_daily_mosaic <- function(date_string) {
  processed$
    filterDate(date_string, as.character(as.Date(date_string) + 1))$
    sort("CLOUDY_PIXEL_PERCENTAGE")$
    mosaic()
}

metadata_path <- file.path(output_dir, paste0("target_image_metadata_", export_label, ".csv"))

if (composite_mode == "monthly_median") {
  export_image <- processed$
    median()$
    clip(study_area_ee)

  if (apply_smoothing == "true") {
    export_image <- export_image$
      focal_mean(
        radius = smoothing_radius,
        kernelType = "square",
        units = "pixels"
      )$
      clip(study_area_ee)
  }

  image_count <- processed$size()$getInfo()

  if (is.null(image_count) || image_count == 0) {
    stop("No Sentinel-2 images found in the requested monthly window within the study area.")
  }

  tibble(
    export_label = export_label,
    start_date = as.character(start_date),
    end_date = as.character(end_date),
    composite_method = "monthly_median",
    n_images = image_count,
    smoothing = apply_smoothing,
    smoothing_radius = smoothing_radius,
    shadow_masking = mask_shadows,
    cloud_threshold = cloud_threshold,
    export_scale = export_scale
  ) %>%
    write_csv(metadata_path)

  message("Composite label: ", export_label)
  message("Composite method: monthly median")
  message("Images used: ", image_count)
  message("Smoothing: ", apply_smoothing)
  message("Smoothing radius: ", smoothing_radius)
  message("Shadow masking: ", mask_shadows)
} else {
  primary_collection <- processed$
    filterDate(target_date, as.character(as.Date(target_date) + 1))

  fill_collection <- processed$
    filterDate(fill_date, as.character(as.Date(fill_date) + 1))

  primary_count <- primary_collection$size()$getInfo()
  fill_count <- fill_collection$size()$getInfo()

  if (is.null(primary_count) || primary_count == 0) {
    stop("No Sentinel-2 image found for target_date within the study area.")
  }

  if (is.null(fill_count) || fill_count == 0) {
    stop("No Sentinel-2 image found for fill_date within the study area.")
  }

  primary_image <- make_daily_mosaic(target_date)
  fill_image <- make_daily_mosaic(fill_date)

  export_image <- primary_image$
    unmask(fill_image)$
    clip(study_area_ee)

  if (apply_smoothing == "true") {
    export_image <- export_image$
      focal_mean(
        radius = smoothing_radius,
        kernelType = "square",
        units = "pixels"
      )$
      clip(study_area_ee)
  }

  primary_id <- tryCatch(
    ee$String(primary_collection$sort("CLOUDY_PIXEL_PERCENTAGE")$first()$get("system:index"))$getInfo(),
    error = function(e) NULL
  )

  fill_id <- tryCatch(
    ee$String(fill_collection$sort("CLOUDY_PIXEL_PERCENTAGE")$first()$get("system:index"))$getInfo(),
    error = function(e) NULL
  )

  tibble(
    export_label = export_label,
    start_date = as.character(start_date),
    end_date = as.character(end_date),
    composite_method = "single_day",
    primary_date = target_date,
    fill_date = fill_date,
    primary_image_id = primary_id,
    fill_image_id = fill_id,
    primary_count = primary_count,
    fill_count = fill_count,
    smoothing = apply_smoothing,
    smoothing_radius = smoothing_radius,
    shadow_masking = mask_shadows,
    cloud_threshold = cloud_threshold,
    export_scale = export_scale
  ) %>%
    write_csv(metadata_path)

  message("Primary date: ", target_date)
  message("Fill date: ", fill_date)
  message("Primary images: ", primary_count)
  message("Fill images: ", fill_count)
  message("Smoothing: ", apply_smoothing)
  message("Smoothing radius: ", smoothing_radius)
  message("Shadow masking: ", mask_shadows)
}

message("Drive folder: ", drive_folder)
message("Export bands: ", paste(export_bands, collapse = ", "))

# geotiff exports
walk(
  export_bands,
  function(band_name) {
    out_path <- file.path(
      output_dir,
      paste0(band_name, "_", export_label, "_mosaic.tif")
    )
    template_path <- file.path(
      paths$raster_dir,
      paste0(band_name, "_", export_label, "_mosaic.tif")
    )
    download_path <- tempfile(
      pattern = paste0(band_name, "_", export_label, "_raw_"),
      tmpdir = output_dir,
      fileext = ".tif"
    )

    ee_as_rast(
      image = export_image$select(band_name),
      dsn = download_path,
      via = "drive",
      container = drive_folder,
      scale = export_scale,
      maxPixels = 1e13,
      timePrefix = FALSE,
      quiet = FALSE
    )

    if (file.exists(template_path)) {
      template_raster <- rast(template_path)
      downloaded_raster <- rast(download_path)

      if (!compareGeom(template_raster, downloaded_raster, stopOnError = FALSE)) {
        downloaded_raster <- project(downloaded_raster, template_raster, method = "bilinear")
      }

      writeRaster(downloaded_raster, out_path, overwrite = TRUE)
    } else {
      file.copy(download_path, out_path, overwrite = TRUE)
    }

    unlink(download_path)
  }
)
