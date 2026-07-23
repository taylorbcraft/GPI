# Load and align Sentinel-2 predictor rasters.
# The helper names each layer, projects or resamples mismatched grids to the
# first predictor, constrains S2REP to its valid range, and returns one stack.

load_predictor_stack <- function(band_names, image_date, resample_method = "near") {
  raster_list <- vector("list", length(band_names))
  names(raster_list) <- band_names
  template <- NULL

  for (band_name in band_names) {
    raster_layer <- rast(path_raster(band_name, image_date))
    names(raster_layer) <- band_name

    # align geometry
    if (is.null(template)) {
      template <- raster_layer
    } else if (!compareGeom(template, raster_layer, stopOnError = FALSE)) {
      if (!same.crs(template, raster_layer)) {
        raster_layer <- project(raster_layer, template, method = resample_method)
      } else {
        raster_layer <- resample(raster_layer, template, method = resample_method)
      }
    }

    # constrain S2REP
    if (band_name == "s2rep") {
      raster_layer <- clamp(raster_layer, lower = 600, upper = 850, values = TRUE)
    }

    raster_list[[band_name]] <- raster_layer
  }

  rast(raster_list)
}
