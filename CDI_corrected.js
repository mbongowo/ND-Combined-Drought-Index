var table = ee.FeatureCollection("projects/ee-geodrought/assets/nd_shape"),
    maxLST = ee.Image("projects/ee-geodrought/assets/maxLST"),
    maxNDVI = ee.Image("projects/ee-geodrought/assets/maxNDVI"),
    minLST = ee.Image("projects/ee-geodrought/assets/minLST"),
    minNDVI = ee.Image("projects/ee-geodrought/assets/minNDVI");
var shp = table;

// Center the map. 
Map.centerObject(shp, 7); 

// Add the ND boundary to the map. 
Map.addLayer(shp, {}, 'North Dakota');

// Set start and end years. 
var startYear = 2023; 
var endYear = 2023; 

// Create two date objects for start and end years. 
var startDate = ee.Date.fromYMD(startYear, 1, 1); 
//print(startDate)
var endDate = ee.Date.fromYMD(endYear + 1, 1, 1);
//print(endDate)

// Make a list with years. 
var years = ee.List.sequence(startYear, endYear); 
//print(years.size(), 'years');

// Make a list with months. 
var months = ee.List.sequence(6, 8);
//print(months.size(), 'months');


var l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(shp)
                          .filterDate('2000-06-01','2011-09-01')
                          .select(['SR_B3', 'SR_B4', 'SR_B5', 'ST_B6', 'QA_PIXEL'])
                          .filter(ee.Filter.calendarRange(6, 8,'month'));
                          
// set  up the cloud masking function
function maskL5srClouds(image) {
  // Collection 2 QA_PIXEL: bits 1, 3, 4 and 5 are dilated cloud, cloud, cloud shadow and snow.
  var dilatedCloudBitMask = (1 << 1);
  var cloudsBitMask = (1 << 3);
  var cloudShadowBitMask = (1 << 4);
  var snowBitMask = (1 << 5);
      // Get the pixel QA band.
  var qa = image.select('QA_PIXEL');
      // All four flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(dilatedCloudBitMask).eq(0)
                    .and(qa.bitwiseAnd(cloudsBitMask).eq(0))
                    .and(qa.bitwiseAnd(cloudShadowBitMask).eq(0))
                    .and(qa.bitwiseAnd(snowBitMask).eq(0));
  return image.updateMask(mask);
}
//SCALE
function L5_applyScaleFactors(image) {
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2).clip(shp);
  var thermalBands = image.select('ST_B6').multiply(0.00341802).add(149.0).clip(shp);
  return image.addBands(opticalBands, null, true)
              .addBands(thermalBands, null, true);
}


var L5= l5.map(L5_applyScaleFactors)
          .map(maskL5srClouds)
          .select(['SR_B3', 'SR_B4', 'SR_B5', 'ST_B6']);



// Define the image collection.
//LANDSAT 7
var L7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(shp)
                      .select(['SR_B3', 'SR_B4', 'SR_B5', 'ST_B6', 'QA_PIXEL'])
                      .filterDate('2012-06-01', '2012-09-01')
                      .filter(ee.Filter.calendarRange(6, 8,'month'));
                      
// set  up the cloud masking function
function maskL7srClouds(image) {
  // Collection 2 QA_PIXEL: bits 1, 3, 4 and 5 are dilated cloud, cloud, cloud shadow and snow.
  var dilatedCloudBitMask = (1 << 1);
  var cloudsBitMask = (1 << 3);
  var cloudShadowBitMask = (1 << 4);
  var snowBitMask = (1 << 5);
      // Get the pixel QA band.
  var qa = image.select('QA_PIXEL');
      // All four flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(dilatedCloudBitMask).eq(0)
                    .and(qa.bitwiseAnd(cloudsBitMask).eq(0))
                    .and(qa.bitwiseAnd(cloudShadowBitMask).eq(0))
                    .and(qa.bitwiseAnd(snowBitMask).eq(0));
  return image.updateMask(mask);
}
//SCALE
function L7_applyScaleFactors(image) {
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2).clip(shp);
  var thermalBands = image.select('ST_B6').multiply(0.00341802).add(149.0).clip(shp);
  return image.addBands(opticalBands, null, true)
              .addBands(thermalBands, null, true);
}

  L7 = L7.map(L7_applyScaleFactors)
                .map(maskL7srClouds)
                .select(['SR_B3', 'SR_B4', 'SR_B5', 'ST_B6']);


//Import and filter the Landsat 8 surface reflectance dataset.
var L8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(shp)
                      .select(['SR_B4', 'SR_B5', 'SR_B6', 'ST_B10', 'QA_PIXEL'])
                      .filter(ee.Filter.calendarRange(6, 8,'month'));

// set  up the cloud masking function
function maskL8srClouds(image) {
  // Collection 2 QA_PIXEL: bits 1, 3, 4 and 5 are dilated cloud, cloud, cloud shadow and snow.
  var dilatedCloudBitMask = (1 << 1);
  var cloudsBitMask = (1 << 3);
  var cloudShadowBitMask = (1 << 4);
  var snowBitMask = (1 << 5);
      // Get the pixel QA band.
  var qa = image.select('QA_PIXEL');
      // All four flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(dilatedCloudBitMask).eq(0)
                    .and(qa.bitwiseAnd(cloudsBitMask).eq(0))
                    .and(qa.bitwiseAnd(cloudShadowBitMask).eq(0))
                    .and(qa.bitwiseAnd(snowBitMask).eq(0));
  return image.updateMask(mask);
}

//Apply scaling factors
function applyScaleFactors(image) {
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2).clip(shp);
  var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0).clip(shp);
  return image.addBands(opticalBands, null, true)
              .addBands(thermalBands, null, true);
}

L8 = L8.map(applyScaleFactors)
          .map(maskL8srClouds)
          .map(function(L8){
            return L8.rename(['SR_B3', 'SR_B4', 'SR_B5', 'ST_B6', 'QA_PIXEL']);
          })
          .select(['SR_B3', 'SR_B4', 'SR_B5', 'ST_B6']); 


var L875 = ee.ImageCollection(L5.merge(L7.merge(L8)));


// The fourteen years analyzed in the paper (Section 2.3.5).
var cdiYears = [2000, 2001, 2002, 2003, 2005, 2006, 2007, 2009, 2011, 2012, 2017, 2018, 2021, 2023];

// Normalization limits (Equation 9) and PCA weights (Section 2.3.6).
// Limits are the 2nd and 98th percentiles of a pooled sample of 5000 cells per
// analyzed year (69,769 cells, 300 m); the weights are the absolute loadings
// of the normalized layers on the first principal component of their
// correlation matrix (58% of the variance). Both come from cdi_corrected_run.py,
// which draws the sample one year at a time: pooling all fourteen years in one
// interactive request fails with "Too many concurrent aggregations".
var limits = {VHI: [0.22868, 0.69866], MSI: [0.45637, 1.42798], WB: [-33.328, 42.411]};
var alpha = 0.9257; // VHI
var beta = 0.9262;  // MSI
var gamma = 0.1300; // WB


// Monthly MSI and VHI for one year, then the mean over the months.
function summerIndices(y) {
  var monthly = ee.ImageCollection.fromImages(
    months.map(function(m) {
      var c = L875.filter(ee.Filter.calendarRange(y, y, 'year'))
                  .filter(ee.Filter.calendarRange(m, m, 'month'));
      return c.mean().set('n', c.size());
    })
  ).filter(ee.Filter.gt('n', 0));   // skip a month with no clear scene

  var indices = monthly.map(function(img) {
    //NDVI
    var NDVI = img.normalizedDifference(['SR_B4', 'SR_B3']);
    //LST
    var LST = img.select('ST_B6');
    //MSI
    var MSI = img.select('SR_B5').divide(img.select('SR_B4')).rename('MSI');
    //VCI
    var VCI = (NDVI.subtract(minNDVI)).divide(maxNDVI.subtract(minNDVI));
    //TCI
    var TCI = (maxLST.subtract(LST)).divide(maxLST.subtract(minLST));
    //VHI (continuous, not classed)
    var VHI = (VCI.multiply(0.5)).add(TCI.multiply(0.5)).rename('VHI');
    return MSI.addBands(VHI);
  });

  return indices.mean().clip(shp);
}


// Import PRISM Precipitation dataset
var PRISM = ee.ImageCollection('OREGONSTATE/PRISM/ANm')
                  .filterDate('2000-01-01', '2024-01-01')
                  .select('ppt');

// Import TerraClimate actual evapotranspiration
var ET = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')
                  .filterDate('2000-01-01', '2024-01-01')
                  .select('aet');

// Monthly water balance (P - ET) for one year, then the mean over the months.
function summerWB(y) {
  var waterBalance = ee.ImageCollection.fromImages(
    months.map(function(m) {
      var P = PRISM.filter(ee.Filter.calendarRange(y, y, 'year'))
        .filter(ee.Filter.calendarRange(m, m, 'month'))
        .sum();
      var E = ET.filter(ee.Filter.calendarRange(y, y, 'year'))
        .filter(ee.Filter.calendarRange(m, m, 'month'))
        .sum()
        .multiply(0.1);
      return P.subtract(E).rename('WB');
    })
  );
  return waterBalance.mean().clip(shp);
}

// VHI, MSI and WB for one year
function layers(y) {
  return summerIndices(y).addBands(summerWB(y));
}

// Equation 9: rescale to [0, 1] with the pooled limits
function normalize(img, lim) {
  return img.subtract(lim[0]).divide(lim[1] - lim[0]).clamp(0, 1);
}

function normLayers(y) {
  var L = layers(y);
  return normalize(L.select('VHI'), limits.VHI).rename('VHI')
    .addBands(normalize(L.select('MSI'), limits.MSI).rename('MSI'))
    .addBands(normalize(L.select('WB'), limits.WB).rename('WB'));
}

// Equation 10: CDI = [alpha(1 - VHI) + beta MSI + gamma(1 - WB)] / 3
function cdiYear(y) {
  var N = normLayers(y);
  return ee.Image(1).subtract(N.select('VHI')).multiply(alpha)
    .add(N.select('MSI').multiply(beta))
    .add(ee.Image(1).subtract(N.select('WB')).multiply(gamma))
    .divide(3)
    .rename('CDI')
    .set('year', y);
}


var cdi = cdiYear(startYear);
Map.addLayer(cdi, {min: 0, max: (alpha + beta + gamma) / 3, palette: ['purple', 'darkblue', 'blue', 'yellow', 'orange', 'red']}, 'Composite Drought Index');

// Masked cells (outside the state, or no clear observation) are written as -9999 and
// flagged as NoData, so they are not read as CDI = 0 (no drought stress).
cdiYears.forEach(function(y) {
  Export.image.toDrive({
    image: cdiYear(y).unmask(-9999),  description: 'CDI_' + y,  folder: 'MSI',  region: shp,  scale: 30,
    crs: 'EPSG:2266',
    maxPixels: 1e13,
    formatOptions: {
      cloudOptimized: true,
      noData: -9999
    }
  });
});
