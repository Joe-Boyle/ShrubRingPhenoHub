// Joe's manuscript: Quick script to extact ERA5 data for QHI
var qhi = ee.Feature(ee.Geometry.Point([-138.90, 69.57]), {name: "qhi"});

// Import and subset ERA5 ImageCollection
var era5 = ee.ImageCollection('ECMWF/ERA5/MONTHLY')
                   .select('total_precipitation')
                   .filter(ee.Filter.date('1990-01-01', '2016-01-01'));

// Extract data for qhi using a reducer
var ERA5_export = era5.map(function(image){ // for each image in the ERA5 collection
  var precip = image.reduceRegion(ee.Reducer.first(), // Reducer first extracts the value of the first pixel overlapping the geometry
                                  qhi.geometry(), // geometry with qhi coordinates
                                  1) // resolution at which to reduce (good practice to set it, here for this purpose set to 1 m)
                    .set('year', image.get('year')) // copy year from image
                    .set('month', image.get('month')); // copy month from image
  return(ee.Feature(qhi.geometry(), precip)); // return as feature as this is required for map
}); 


// Export FC as table to drive. 
Export.table.toDrive({
  collection: ERA5_export,
  description: 'ERA_export_joe_QHI',
  fileNamePrefix: 'ERA_export_joe_QHI',
  fileFormat: 'csv'
});