/*
  Moisture Stress Index (MSI), North Dakota, Google Earth Engine
  Supplementary code: Mbuh, Caparo and Alu, Mapping drought patterns in North Dakota, USA,
  with an experimental Landsat-based Combined Drought Index

  MSI = SWIR1/NIR (SR_B6/SR_B5), Landsat 8 Collection 2 Level 2 surface reflectance.
  QA_PIXEL bits 3 and 5 are masked (cloud and snow in Collection 2).
  Monthly means for May to September; annual MSI is the mean of the monthly means.
  Display stretch 0.25 to 1.0 (dark blue to red), legend classes Normal to Extreme.

  Changed from our working copy so it runs outside our account:
    ND boundary from TIGER/2018/States in place of our uploaded asset (table)
    block comment markers restored
    var added to MSI2000
*/
var shp = ee.FeatureCollection('TIGER/2018/States').filter(ee.Filter.eq('NAME', 'North Dakota'));
// Center the map.
Map.centerObject(shp, 6.5);
// Add the ND boundary to the map.
Map.addLayer(shp, {}, 'North Dakota');
// Set start and end years.
var startYear = 2018;
var endYear = 2023;
// Create two date objects for start and end years.
var startDate = ee.Date.fromYMD(startYear, 1, 1);
print(startDate)
var endDate = ee.Date.fromYMD(endYear + 1, 1, 1);
print(endDate)
// Make a list with years.
var years = ee.List.sequence(startYear, endYear);
//print(years.size(), 'years');
// Make a list with months.
var months = ee.List.sequence(5, 9);
//print(months.size(), 'months');
//Import and filter the Landsat 8 surface reflectance dataset.
var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(shp)
.filter(ee.Filter.date(startDate, endDate));
// set up the cloud masking function
function maskL8srClouds(image) {
// Bits 3 and 5 are cloud shadow and cloud, respectively.
var cloudShadowBitMask = (1 << 3);
var cloudsBitMask = (1 << 5);
// Get the pixel QA band.
var qa = image.select('QA_PIXEL');
// Both flags should be set to zero, indicating clear conditions.
var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
.and(qa.bitwiseAnd(cloudsBitMask).eq(0));
return image.updateMask(mask);
}
//Apply scaling factors
function applyScaleFactors(image) {
var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2).clip(shp);
var thermalBands = image.select('ST_B.').multiply(0.00341802).add(149.0).clip(shp);
return image.addBands(opticalBands, null, true)
.addBands(thermalBands, null, true);
}
var lsat8 = l8.map(applyScaleFactors)
.map(maskL8srClouds);
// We use a function to calculate the Moisture Stress Index.
// We map over the L8 image collection and select the NIR and SWIR bands
// We set the timestamp and return the MSI.
var MSI = lsat8.map(function(image) {
var nirband = image.select('SR_B5');
var swirband = image.select('SR_B6');
var msi = swirband.divide(nirband).rename('MSI')
.set('system:time_start', image.get(
'system:time_start'));
return msi;
});
// We apply a nested loop where we first map over the relevant years and then map over the relevant months.
//The function returns an image with MSI for each month.
// A flatten is applied to convert a feature collection of features into a single feature collection.
var monthlyMSI = ee.ImageCollection.fromImages(
years.map(function(y) {
return months.map(function(m) {
var msi = MSI.filter(ee.Filter
.calendarRange(y, y, 'year'))
.filter(ee.Filter.calendarRange(m, m, 'month'))
.mean();
return msi.set('year', y)
.set('month', m)
.set('system:time_start', ee.Date
.fromYMD(y, m, 1));
});
}).flatten()
);
// Add the layer with monthly mean. Clip for the ND State boundary.
var moisture_stress_index = MSI.select('msi');
var msiVis = {
min: 0.25,
max: 1,
palette: ['darkblue', 'blue', 'yellow', 'orange', 'red'],
};
Map.addLayer(monthlyMSI.mean().clip(shp),
msiVis,
'6 years Mean MSI');
var MSI_2018 = monthlyMSI.filterDate('2018-05-01', '2019-05-01').mean().clip(shp)
.reproject({'crs': 'EPSG:32613', 'crsTransform': null, 'scale': 30});
var MSI_2019 = monthlyMSI.filterDate('2019-05-01', '2020-05-01').mean().clip(shp)
.reproject({'crs': 'EPSG:32613', 'crsTransform': null, 'scale': 30});
var MSI_2020 = monthlyMSI.filterDate('2020-05-01', '2021-05-01').mean().clip(shp)
.reproject({'crs': 'EPSG:32613', 'crsTransform': null, 'scale': 30});
print (MSI_2020);
var MSI_2021 = monthlyMSI.filterDate('2021-05-01', '2022-05-01').mean().clip(shp)
.reproject({'crs': 'EPSG:32613', 'crsTransform': null, 'scale': 30});
var MSI_2022 = monthlyMSI.filterDate('2022-05-01', '2023-05-01').mean().clip(shp)
.reproject({'crs': 'EPSG:32613', 'crsTransform': null, 'scale': 30});
var MSI_2023 = monthlyMSI.filterDate('2023-05-01', '2024-05-01').mean().clip(shp)
.reproject({'crs': 'EPSG:32613', 'crsTransform': null, 'scale': 30});
//print(MSI_2018)
//FIND MSI 2018 SCALE
var prop= MSI_2018.select('MSI').projection()//.nominalScale();
//print (prop);
var projection = MSI_2018.select('MSI').projection().getInfo();
//print(projection);
//Map.addLayer(MSI_2018, msiVis, 'MSI2018');
/*Map.addLayer(monthlyMSI.filterDate('2019-05-01', '2020-05-01').mean().clip(shp),
msiVis,
'MSI2019');*/
Map.addLayer(monthlyMSI.filterDate('2020-05-01', '2021-05-01').mean().clip(shp),
msiVis,
'MSI2020');
/*Map.addLayer(monthlyMSI.filterDate('2021-05-01', '2022-05-01').mean().clip(shp),
msiVis,
'MSI2021');
Map.addLayer(monthlyMSI.filterDate('2022-05-01', '2023-05-01').mean().clip(shp),
msiVis,
'MSI2022');*/
var MSI2000 = monthlyMSI.filterDate('2023-05-01', '2024-05-01').mean().clip(shp);
//MSI Legend
///////////////////////////////////////////////////////////////////////////////
var legend = ui.Panel({
style: {
position: 'bottom-left',
padding: '8px 15px'
}
});
// Create and add the legend title.
var legendTitle = ui.Label({
value: 'MSI',
style: {
fontWeight: 'bold',
fontSize: '18px',
margin: '0 0 4px 0',
padding: '0'
}
});
legend.add(legendTitle);
// Creates and styles 1 row of the legend.
var makeRow = function(color, name) {
// Create the label that is actually the colored box.
var colorBox = ui.Label({
style: {
backgroundColor: '#' + color,
// Use padding to give the box height and width.
padding: '8px',
margin: '0 0 4px 0'
}
});
// Create the label filled with the description text.
var description = ui.Label({
value: name,
style: {margin: '0 0 4px 6px'}
});
return ui.Panel({
widgets: [colorBox, description],
layout: ui.Panel.Layout.Flow('horizontal')
});
};
// Legend Rows
legend.add(makeRow('00008b', 'Normal'));//1
legend.add(makeRow('0000ff', 'Mild'));//2
legend.add(makeRow('ffff00', 'Moderate'));//3
legend.add(makeRow('ffa500', 'Severe'));//4
legend.add(makeRow('ff0000', 'Extreme'));//5
// Add Legend to Map
Map.add(legend);
// Set the title and axis labels for the chart.
var title = {
title: 'Monthly MSI',
hAxis: {
title: 'Time'
},
vAxis: {
title: 'MSI'
},
colors: ['darkgreen']
};
// Plot the chart using the ND boundary.
var chartMonthly = ui.Chart.image.seriesByRegion({
imageCollection: monthlyMSI,
regions: shp.geometry(),
reducer: ee.Reducer.mean(),
band: 'MSI',
scale: 30,
xProperty: 'system:time_start'
}).setSeriesNames(['MSI'])
.setOptions(title)
.setChartType('ColumnChart');
// Print the chart.
print(chartMonthly);
/*Export.image.toDrive({
image: MSI_2023, description: 'MSI_2023', folder: 'MSI', region: shp, scale: 30,
crs: projection.crs,
crsTransform: projection.transform,
//crsTransform: projection.transform,
maxPixels: 1e13,
// formatOptions: {
// cloudOptimized: true
// }
});*/
/*var batch = require('users/fitoprincipe/geetools:batch')
//var p = ee.Geometry.Point([-72, -42])
var col = ee.ImageCollection('COPERNICUS/S2')
.filterBounds(shp)

var asset = 'ic'
var options = {
name: 'msi_{system_date}',
scale: 30,
maxPixels: 1e13,
region: shp
}
batch.Download.ImageCollection.toAsset(col, asset, options)*/
