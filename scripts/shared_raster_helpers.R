# Load Sentinel-2 predictor rasters.
# The helper names each layer, constrains S2REP to its valid range, and returns
# one stack.

load_predictor_stack <- function(band_names, image_date) {
  raster_list <- vector("list", length(band_names))
  names(raster_list) <- band_names

  for (band_name in band_names) {
    raster_layer <- rast(path_raster(band_name, image_date))
    names(raster_layer) <- band_name

    # constrain S2REP
    if (band_name == "s2rep") {
      raster_layer <- clamp(raster_layer, lower = 600, upper = 850, values = TRUE)
    }

    raster_list[[band_name]] <- raster_layer
  }

  rast(raster_list)
}
