// Record and rebuild the Sentinel-2 predictor rasters used by the R pipeline.
// Run this in Google Earth Engine only when data/processed/rasters needs to be
// regenerated. The local workflow expects one GeoTIFF per predictor band named
// <band>_<date>_mosaic.tif.

// Settings
var studyArea = geometry;

var year = 2025;
var startDate = '2025-04-01';
var endDate = '2025-04-30';
var targetDate = '2025-04-11';
var cloudThreshold = 65;
var exportScale = 20;
var exportFolder = 'gpi_2025';
var exportPrefixDate = '2025-04-11';

// Sentinel-2 surface reflectance and cloud probability collections

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

// Cloud mask, reflectance scaling, and spectral indices

function maskClouds(img) {
  var cloudProb = ee.Image(img.get('cloud_mask')).select('probability');
  var cloudMask = cloudProb.lt(cloudThreshold);

  return img
    .updateMask(cloudMask)
    .divide(10000)
    .copyProperties(img, img.propertyNames());
}

function addIndices(img) {
  var s2rep = img.expression(
    '705 + 35 * ((((R + RE3) / 2) - RE1) / (RE2 - RE1))',
    {
      R: img.select('B4'),
      RE1: img.select('B5'),
      RE2: img.select('B6'),
      RE3: img.select('B7')
    }
  ).rename('s2rep');

  var ndvi = img.normalizedDifference(['B8', 'B4']).rename('ndvi');
  var ndwi = img.normalizedDifference(['B3', 'B8']).rename('ndwi');
  var savi = img.expression(
    '((NIR - RED) / (NIR + RED + L)) * (1 + L)',
    {
      NIR: img.select('B8'),
      RED: img.select('B4'),
      L: 0.5
    }
  ).rename('savi');
  var evi = img.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))',
    {
      NIR: img.select('B8'),
      RED: img.select('B4'),
      BLUE: img.select('B2')
    }
  ).rename('evi');
  var msi = img.expression(
    'SWIR1 / NIR',
    {
      SWIR1: img.select('B11'),
      NIR: img.select('B8')
    }
  ).rename('msi');
  var ndmi = img.normalizedDifference(['B8', 'B11']).rename('ndmi');
  var mndwi = img.normalizedDifference(['B3', 'B11']).rename('mndwi');

  return img.addBands([s2rep, ndvi, ndwi, savi, evi, msi, ndmi, mndwi]);
}

// Build the target image

var processed = joinedCollection
  .map(maskClouds)
  .map(addIndices);

var targetImage = processed
  .filterDate(targetDate, ee.Date(targetDate).advance(1, 'day'))
  .sort('CLOUDY_PIXEL_PERCENTAGE')
  .first();

targetImage = ee.Image(targetImage).clip(studyArea);

// Visual checks in the Earth Engine map

print('target image', targetImage);
print('target image id', targetImage.get('system:index'));
print('target image date', targetImage.date());

Map.centerObject(studyArea, 13);
Map.addLayer(
  targetImage,
  {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.3},
  'rgb'
);
Map.addLayer(
  targetImage.select('s2rep'),
  {min: 710, max: 730, palette: ['yellow', 'green', 'darkgreen']},
  's2rep'
);
Map.addLayer(
  targetImage.select('evi'),
  {min: 0, max: 1, palette: ['white', 'lightgreen', 'darkgreen']},
  'evi'
);

// Export the raster stack expected by config.R

function exportBand(image, bandName, fileName) {
  Export.image.toDrive({
    image: image.select(bandName),
    description: fileName.replace('.tif', ''),
    folder: exportFolder,
    fileNamePrefix: fileName.replace('.tif', ''),
    region: studyArea,
    scale: exportScale,
    maxPixels: 1e13
  });
}

exportBand(targetImage, 's2rep', 's2rep_' + exportPrefixDate + '_mosaic.tif');
exportBand(targetImage, 'ndvi',  'ndvi_'  + exportPrefixDate + '_mosaic.tif');
exportBand(targetImage, 'ndwi',  'ndwi_'  + exportPrefixDate + '_mosaic.tif');
exportBand(targetImage, 'savi',  'savi_'  + exportPrefixDate + '_mosaic.tif');
exportBand(targetImage, 'evi',   'evi_'   + exportPrefixDate + '_mosaic.tif');
exportBand(targetImage, 'msi',   'msi_'   + exportPrefixDate + '_mosaic.tif');
exportBand(targetImage, 'ndmi',  'ndmi_'  + exportPrefixDate + '_mosaic.tif');
exportBand(targetImage, 'mndwi', 'mndwi_' + exportPrefixDate + '_mosaic.tif');
