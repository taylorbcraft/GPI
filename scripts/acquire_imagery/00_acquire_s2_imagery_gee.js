// Acquire Sentinel-2 predictor rasters in the Earth Engine Code Editor.
// Update the settings below, then run the script and start the tasks in Exports.

// settings
var startDate = '2025-04-01';
var endDate = '2025-05-01'; // exclusive
var exportLabel = '2025-04-median';
var applySmoothing = true;
var smoothingRadius = 3;
var maskShadows = false;
var cloudThreshold = 65;
var exportScale = 20;
var driveFolder = 'GEE';
var exportBands = [
  's2rep',
  'ndvi',
  'ndwi',
  'savi',
  'evi',
  'msi',
  'ndmi',
  'mndwi'
];

// field geometry bounding box in EPSG:4326
var studyArea = ee.Geometry.Rectangle(
  [5.359333415987, 52.848757420800, 5.603293926403, 53.055335338008],
  'EPSG:4326',
  false
);

// source collections
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(studyArea)
  .filterDate(startDate, endDate);

var s2Clouds = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
  .filterBounds(studyArea)
  .filterDate(startDate, endDate);

var joined = ee.Join.saveFirst('cloud_mask').apply({
  primary: s2,
  secondary: s2Clouds,
  condition: ee.Filter.equals({
    leftField: 'system:index',
    rightField: 'system:index'
  })
});

var joinedCollection = ee.ImageCollection(joined);

// cloud mask
var maskClouds = function(image) {
  image = ee.Image(image);
  var cloudProbability = ee.Image(image.get('cloud_mask'))
    .select('probability');
  var combinedMask = cloudProbability.lt(cloudThreshold);

  if (maskShadows) {
    combinedMask = combinedMask.and(image.select('SCL').neq(3));
  }

  return image
    .updateMask(combinedMask)
    .divide(10000)
    .copyProperties(image, image.propertyNames());
};

// predictor bands
var addIndices = function(image) {
  var s2rep = image.expression(
    '705 + 35 * ((((R + RE3) / 2) - RE1) / (RE2 - RE1))',
    {
      R: image.select('B4'),
      RE1: image.select('B5'),
      RE2: image.select('B6'),
      RE3: image.select('B7')
    }
  ).rename('s2rep');

  var ndvi = image.normalizedDifference(['B8', 'B4']).rename('ndvi');
  var ndwi = image.normalizedDifference(['B3', 'B8']).rename('ndwi');

  var savi = image.expression(
    '((NIR - RED) / (NIR + RED + L)) * (1 + L)',
    {
      NIR: image.select('B8'),
      RED: image.select('B4'),
      L: 0.5
    }
  ).rename('savi');

  var evi = image.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))',
    {
      NIR: image.select('B8'),
      RED: image.select('B4'),
      BLUE: image.select('B2')
    }
  ).rename('evi');

  var msi = image.expression('SWIR1 / NIR', {
    SWIR1: image.select('B11'),
    NIR: image.select('B8')
  }).rename('msi');

  var ndmi = image.normalizedDifference(['B8', 'B11']).rename('ndmi');
  var mndwi = image.normalizedDifference(['B3', 'B11']).rename('mndwi');

  return image
    .addBands(s2rep)
    .addBands(ndvi)
    .addBands(ndwi)
    .addBands(savi)
    .addBands(evi)
    .addBands(msi)
    .addBands(ndmi)
    .addBands(mndwi);
};

var processed = joinedCollection
  .map(maskClouds)
  .map(addIndices);

// build composite
var exportImage = processed.median().clip(studyArea);

// optional spatial smoothing
if (applySmoothing) {
  exportImage = exportImage.focalMean({
    radius: smoothingRadius,
    kernelType: 'square',
    units: 'pixels'
  }).clip(studyArea);
}

print('Composite label', exportLabel);
print('Composite method', 'monthly median');
print('Images used', processed.size());
print('Smoothing', applySmoothing);
print('Smoothing radius', smoothingRadius);
print('Shadow masking', maskShadows);
print('Export bands', exportBands);

Map.centerObject(studyArea, 10);
Map.addLayer(
  exportImage.select('s2rep'),
  {min: 700, max: 750, palette: ['440154', '21918c', 'fde725']},
  'S2REP ' + exportLabel
);

// geotiff exports
exportBands.forEach(function(bandName) {
  var exportName = bandName + '_' + exportLabel + '_mosaic';

  Export.image.toDrive({
    image: exportImage.select(bandName),
    description: exportName,
    folder: driveFolder,
    fileNamePrefix: exportName,
    region: studyArea,
    scale: exportScale,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  });
});

// export run metadata
var metadataProperties = {
  export_label: exportLabel,
  start_date: startDate,
  end_date: ee.Date(endDate).advance(-1, 'day').format('YYYY-MM-dd'),
  composite_method: 'monthly_median',
  n_images: processed.size(),
  smoothing: applySmoothing,
  smoothing_radius: smoothingRadius,
  shadow_masking: maskShadows,
  cloud_threshold: cloudThreshold,
  export_scale: exportScale
};

var metadata = ee.FeatureCollection([
  ee.Feature(null, metadataProperties)
]);

Export.table.toDrive({
  collection: metadata,
  description: 's2_composite_metadata_' + exportLabel,
  folder: driveFolder,
  fileNamePrefix: 's2_composite_metadata_' + exportLabel,
  fileFormat: 'CSV'
});
