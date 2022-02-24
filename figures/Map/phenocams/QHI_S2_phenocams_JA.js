// QHI Phenocam Sentinel Time-Series extractions for Joe
// Jakob Assmann March 2021 j.assmann@bio.au.dk

// Set phenocam locations
var phenocams = ee.FeatureCollection([
  ee.Feature(ee.Geometry.Point([-138.90460, 69.57560]), {name: 'Phenocam1'}),
  ee.Feature(ee.Geometry.Point([-138.90525, 69.57570]), {name: 'Phenocam2'}),
  ee.Feature(ee.Geometry.Point([-138.91277, 69.57781]), {name: 'Phenocam3'}),
  ee.Feature(ee.Geometry.Point([-138.89453, 69.57515]), {name: 'Phenocam4'}),
  ee.Feature(ee.Geometry.Point([-138.86711, 69.57651]), {name: 'Phenocam5'}),
  ee.Feature(ee.Geometry.Point([-138.86305, 69.57483]), {name: 'Phenocam6'})]);

// Add coordinates as properties for convenience post export
phenocams = phenocams.map(function(feature){
  var long = feature.geometry().coordinates().get(0);
  var lat =feature.geometry().coordinates().get(1);
  return(feature.set('long', long)
                .set('lat', lat));
});

// Check location data
print(phenocams, 'Phenocam Locations');
Map.addLayer(phenocams, {name :'Phenocams', color: "red"});
Map.centerObject(phenocams);

// Filter sentinel data (Level 1, 2 and cloud probabillity)
var s2l1_phenocam = s2l1
  .filterBounds(phenocams)
  .filterDate("2016-01-01", "2019-03-01");
print(s2l1_phenocam.size(), 'Number of S2 Level 1 Images');
var s2l2_phenocam = s2l2
  .filterBounds(phenocams)
  .filterDate("2016-01-01", "2019-03-01");
print(s2l2_phenocam.size(), 'Number of S2 Level 2 Images');
var s2cloudless_phenocam = s2cloudless
  .filterBounds(phenocams)
  .filterDate("2016-01-01", "2019-03-01");
print(s2cloudless_phenocam.size(), 'Number of S2 Cloud Propability Esimtates');

// Only 10 Level 2 images, we'll just focus on Level 1

// Notes on cloud probabliltiy 's2cloudless'
// The probabilities 0-100 % are generated with https://github.com/sentinel-hub/sentinel2-cloud-detector
// See also https://medium.com/sentinel-hub/cloud-masks-at-your-service-6e5b2cb2ce8a
// The values stored in the GEE a cloud probabilties for 
// the 10 m x 10 m pixels (higher res bands interpolated to 10 m)

// Define function to add indices
var addIndices = function(s2Image){
  // NDVI and kNDVI at 10 m
  var s2NDVI = s2Image.normalizedDifference(['B8','B4']).rename('NDVI');
  var s2kNDVI  = s2NDVI.pow(2).tanh().rename('kNDVI');

  // NDSI
  // Resample B3 (10 m) to resolution fo B11 (20 m) using bicubic interpolation
  var s2B11 = s2Image.select('B11');
  var s2B3_20m = s2Image.select('B3').resample('bicubic')
                          .reproject(s2B11.projection()).rename('B3_20m');
  var s2B3B11_20m = s2B3_20m.addBands(s2B11);
  var s2NDSI_20m = s2B3B11_20m.normalizedDifference(['B3_20m','B11']).rename('NDSI_20m');

  // NDVI and kNDVI at 20 m
  // Resample bands for direct comparability with NDSI (20 m)
  var s2B4_20m = s2Image.select('B4').resample('bicubic')
                          .reproject(s2B11.projection()).rename('B4_20m');
  var s2B8_20m = s2Image.select('B8').resample('bicubic')
                          .reproject(s2B11.projection()).rename('B8_20m');
  var s2B4B8_20m = s2B4_20m.addBands(s2B8_20m);
  var s2NDVI_20m = s2B4B8_20m.normalizedDifference(['B8_20m','B4_20m']).rename('NDVI_20m');
  var s2kNDVI_20m  = s2NDVI_20m.pow(2).tanh().rename('kNDVI_20m');
  
  // Combine all bands and indices into new image and return
  return(s2Image.addBands(s2B3_20m)
                .addBands(s2B4_20m)
                .addBands(s2B8_20m)
                .addBands(s2NDVI)
                .addBands(s2kNDVI)
                .addBands(s2NDVI_20m)
                .addBands(s2kNDVI_20m)
                .addBands(s2NDSI_20m)
                );
};

// Define function to add cloud propability
// This sections causes three errors that I cannot explain
// but these also do not seem to impact the calculations 
var addCloudProp = function(s2Image){
  var s2Image_systemID = s2Image.get("system:index");
  print(s2Image_systemID);
  var cloudProp = s2cloudless_phenocam
                        .filterMetadata('system:index', 'equals', s2Image_systemID)
                        .first().rename("CLP");
  var cloduProp_20m = cloudProp.select('CLP').resample('bicubic')
                          .reproject(s2Image.select('B11').projection()).rename('CLP_20m')
  return(s2Image.addBands(cloudProp)
                .addBands(cloduProp_20m));
};

// Add indices, cloud propabilities and image stats to s2 image collections
var s2l1_phenocam_final = s2l1_phenocam.map(function(image){
  var year = ee.Image(ee.Number.parse(ee.Date(image.get('system:time_start')).format('y'))).rename('year');
  var doy = ee.Image(ee.Number.parse(ee.Date(image.get('system:time_start')).format('D'))).rename('doi');
  var projection = ee.Image(ee.Number.parse(image.select('B1').projection().crs().slice(5))).rename('EPSG');
  image = image.addBands(year).addBands(doy).addBands(projection);
  image = addIndices(image);
  image = addCloudProp(image);
  image = image.select(['year',
                       'doi',
                       'EPSG',
                       'CLP',
                       'CLP_20m',
                       'B3',
                       'B3_20m',
                       'B4',
                       'B4_20m',
                       'B8',
                       'B8_20m',
                       'B11',
                       'NDVI',
                       'NDVI_20m',
                       'kNDVI',
                       'kNDVI_20m',
                       'NDSI_20m']);
  return(image);
});

// Extract variables for phenocam sites
var phenocam_time_series = s2l1_phenocam_final.map(
  function(image){
    return(image.reduceRegions({
        collection: phenocams,
        reducer: ee.Reducer.first(),
        scale: 10,
        crs: 'EPSG:32607'
      }));
}).flatten();

// Print first objedt to check everythig worked
print(phenocam_time_series.first());

// Export to Drive
Export.table.toDrive(phenocam_time_series,
  'JoesQHIphenocamExport',
  '',
  'S2QHIphenocam',
  'csv');
  
// Print a message on the console
print("Joe is hot! Even Metro says so:");
print("https://metro.co.uk/2016/09/13/ireland-or-aland-university-challenge-contestants-answer-to-island-question-divides-viewers-6125265/");