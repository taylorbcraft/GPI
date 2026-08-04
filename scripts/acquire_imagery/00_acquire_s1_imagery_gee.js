// Build a Sentinel-1 temporal variability raster in the Earth Engine Code Editor.
// Update the settings below, then run the script and start the task in Exports.

// settings
var s1Year = '2025';
var startDate = '2025-03-31';
var endDate = '2025-07-17'; // inclusive
var orbitPass = 'ASCENDING';
var relativeOrbit = 88;
var polarization = 'VV';
var instrumentMode = 'IW';
var boxcarRadius = 3;
var exportScale = 10;
var driveFolder = 'GEE';

// field geometry bounding box in EPSG:4326
var studyArea = ee.Geometry.Rectangle(
  [5.359333415987, 52.848757420800, 5.603293926403, 53.055335338008],
  'EPSG:4326',
  false
);

// filter Sentinel-1 scenes
var s1Collection = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(studyArea)
  .filterDate(startDate, ee.Date(endDate).advance(1, 'day'))
  .filter(ee.Filter.eq('orbitProperties_pass', orbitPass))
  .filter(ee.Filter.inList('relativeOrbitNumber_start', [relativeOrbit]))
  .filter(ee.Filter.listContains(
    'transmitterReceiverPolarisation',
    polarization
  ))
  .filter(ee.Filter.eq('instrumentMode', instrumentMode))
  .select(polarization)
  .sort('system:time_start');

// spatial smoothing in linear power
var boxcar = ee.Kernel.square({
  radius: boxcarRadius,
  units: 'pixels'
});

var s1LinearSmoothed = s1Collection.map(function(image) {
  var linearImage = ee.Image.constant(10).pow(image.divide(10));

  return linearImage
    .convolve(boxcar)
    .rename('linear_smoothed')
    .copyProperties(image, image.propertyNames());
});

// consecutive image pairs
var imageList = s1LinearSmoothed.toList(s1LinearSmoothed.size());
var logRatioList = ee.List.sequence(
  0,
  s1LinearSmoothed.size().subtract(2)
).map(function(index) {
  var i = ee.Number(index);
  var image1 = ee.Image(imageList.get(i));
  var image2 = ee.Image(imageList.get(i.add(1)));

  return image2.log10()
    .subtract(image1.log10())
    .rename('log_ratio')
    .set('system:time_start', image2.get('system:time_start'));
});

var logRatioCollection = ee.ImageCollection.fromImages(logRatioList);

// reduce temporal variation
var temporalSd = logRatioCollection
  .reduce(ee.Reducer.stdDev())
  .clip(studyArea);

var exportName = [
  'S1',
  polarization,
  'LogRatio_StdDev',
  orbitPass.charAt(0) + orbitPass.slice(1).toLowerCase(),
  s1Year
].join('_');

print('Images used', s1Collection.size());
print('Consecutive log ratios', logRatioCollection.size());
print('Orbit pass', orbitPass);
print('Relative orbit', relativeOrbit);

Map.centerObject(studyArea, 10);
Map.addLayer(
  temporalSd,
  {min: 0, max: 0.15, palette: ['ffffff', 'fdae61', 'd7191c']},
  exportName
);

// export raster
Export.image.toDrive({
  image: temporalSd,
  description: exportName,
  folder: driveFolder,
  fileNamePrefix: exportName,
  region: studyArea,
  scale: exportScale,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

// export run metadata
var metadata = ee.FeatureCollection([
  ee.Feature(null, {
    s1_year: s1Year,
    start_date: startDate,
    end_date: endDate,
    orbit_pass: orbitPass,
    relative_orbit: relativeOrbit,
    polarization: polarization,
    instrument_mode: instrumentMode,
    image_count: s1Collection.size(),
    log_ratio_count: logRatioCollection.size(),
    boxcar_radius: boxcarRadius,
    export_scale: exportScale
  })
]);

Export.table.toDrive({
  collection: metadata,
  description: 's1_temporal_sd_metadata_' + s1Year,
  folder: driveFolder,
  fileNamePrefix: 's1_temporal_sd_metadata_' + s1Year,
  fileFormat: 'CSV'
});
