var s2 = ee.ImageCollection("COPERNICUS/S2"),
    county = ee.FeatureCollection("***"), //Administrative Divisions
    fenge_liuan = ee.FeatureCollection("***"); //Segmented Objects

var rebuildImageLib = require('***/defult:lib/rebuildImageLib');
var commonLib = require("***/defult:lib/commonLib");

var pixelScale = 10;
var rawImageStartDate = "2017-3-1";
var rawImageEndDate = "2017-10-31";
var treeStartDate = "2017-3-1";
var treeEndDate = "2017-5-1";

//evi2Tilering timing
var selectStartDate = "2017-7-1";
var selectEndDate = "2017-8-20";

var govAreaCode = 341502;
var selectBandList = ["EVI2", "LSWI", "NDVI", "NDWI", "NDBI"];

//shapefile boulder
var govArea = county.filter(ee.Filter.eq("ADMINCODE", govAreaCode.toString()));
Map.centerObject(govArea, 9);
//masks
function maskCollection(segRawCol) {
  //water mask
  var maxNDWI = segRawCol.select("NDWI").reduce(ee.Reducer.max());
  var waterMask = maxNDWI.lte(0.2);
  Map.addLayer(waterMask.updateMask(waterMask), {palette:"0000ff"}, "noWater", false);
  
  //urban mask
  var minNDBI = segRawCol.select("NDBI").reduce(ee.Reducer.min());
  var urbanMask = minNDBI.lt(-0.2);
  Map.addLayer(urbanMask.updateMask(urbanMask), {palette:"00ff00"}, "noUrban", false);
  
  //ndvi mask
  var minMaxNDVI = segRawCol.select("NDVI").reduce(ee.Reducer.minMax());
  var maxNDVI = minMaxNDVI.select("NDVI_max");
  var minNDVI = minMaxNDVI.select("NDVI_min");
  var ndviMask = maxNDVI.gte(0.57)
                        .and(maxNDVI.subtract(minNDVI).gt(0.5));
  Map.addLayer(ndviMask.updateMask(ndviMask), {palette:"ffff00"}, "greenNDVI", false);
  
  //tree mask
  var segRawCol2 = segRawCol.filterDate(treeStartDate, treeEndDate);
  var treeMinMaxNDVI = segRawCol2.select("NDVI").reduce(ee.Reducer.minMax());
  var treeMaxNDVI = treeMinMaxNDVI.select("NDVI_max");
  var treeMinNDVI = treeMinMaxNDVI.select("NDVI_min");
  var treeNDVIMask = treeMaxNDVI.lte(0.5);
  Map.addLayer(treeNDVIMask.updateMask(treeNDVIMask), {palette:"aabbcc"}, "notTreeNDVI", false);

  var s2Collection = segRawCol.map(function(img){
                                img = img.updateMask(urbanMask);
                                img = img.updateMask(waterMask);
                                img = img.updateMask(ndviMask);
                                img = img.updateMask(treeNDVIMask);
                                return img;
                              });
  return s2Collection;
}

//Indecies Calculation
function getSegmentResult(roi, bandNames, rawCol) {
  var segRawCol = maskCollection(rawCol);
  print("segRawCol", segRawCol);
  var evi2Heading = segRawCol.select("EVI2").reduce(ee.Reducer.max());
  var evi2HeadingImg = segRawCol.map(function(img){
    img = img.updateMask(img.select("EVI2").eq(evi2Heading));
    img = img.addBands(img.select("EVI2").rename("EVI2_heading"));
    img = img.addBands(img.select("LSWI").rename("LSWI_max"));
    return img;
  }).mosaic();
  var lswiMax = evi2HeadingImg.select("LSWI_max");
  
  var segRawCol3 = segRawCol.filterDate(selectStartDate, selectEndDate);
  var evi2Tilering = segRawCol3.select("EVI2").reduce(ee.Reducer.min());
  var evi2TileringImg = segRawCol3.map(function(img){
    img = img.updateMask(img.select("EVI2").eq(evi2Tilering));
    img = img.addBands(img.select("EVI2").rename("EVI2_tilering"));
    return img.addBands(img.select("LSWI").rename("LSWI_min"));
  }).mosaic();
  var lswiMin = evi2TileringImg.select("LSWI_min");
  
  var img1 = ee.Image.constant(0).clip(roi);
  img1 = img1.addBands(evi2Heading.rename("EVI2_heading"));
  img1 = img1.addBands(lswiMax.rename("LSWI_max"));
  img1 = img1.addBands(evi2Tilering.rename("EVI2_tilering"));
  img1 = img1.addBands(lswiMin.rename("LSWI_min"));
  var rcleImage = rebuildImageLib.s2RCLE(img1);
  rcleImage = rcleImage.addBands(ee.Image.constant(1).rename("result"));
  var mask = rcleImage.select("LSWI_min").gt(0.1)
                      .and(rcleImage.select("RCLE").lt(0.6));
  var resultImage = rcleImage.updateMask(mask).select("result");
  return resultImage;
}

// process segmented images
function processS2Data(startDate, endDate, bandNames) {
  var s2RawCol = rebuildImageLib.getS2Collection(startDate, endDate, govArea, bandNames);
  // var roi = fenge_liuan.filter(ee.Filter.and(ee.Filter.gte("index", 37959), ee.Filter.lte("index", 38263)));
  ////////rebuild imageCollection based on the Objects////////
  var roi = fenge_liuan;
  var selectRawCol = s2RawCol.map(function(img){
    var fCol = img.clip(govArea)
                  .select(bandNames)
                  .reduceRegions({
                    reducer: ee.Reducer.mean(),
                    collection: roi,
                    scale: pixelScale
                  });
    
    var img1 = ee.Image();
    img1 = img1.set("system:index", img.get("system:index"));
    img1 = img1.set("system:time_start", img.get("system:time_start"));
    for (var i=0; i<bandNames.length; i++) {
      var bandName = bandNames[i];
      var bandImg = fCol.reduceToImage([bandName], ee.Reducer.mean());
      img1 = img1.addBands(bandImg.rename(bandName));
    }
    img1 = img1.clip(govArea);
    return img1;
  }).select(bandNames);
  ////////////////////////////////////////////
  
  //result
  var segmentResult = getSegmentResult(roi, bandNames, selectRawCol);
  var result = ee.Image.constant(0).clip(roi);
  result = result.where(segmentResult.eq(1), 1);
  print(result);
  Map.addLayer(result.updateMask(result), {palette:"ff0000"}, "result");
    Export.image.toDrive({
      image: result, //Export result
      description: "result-q",  
      fileNamePrefix: "result-q", //Export File Name
      crs: "EPSG:4326", //proj info
      region: govArea, //area
      scale: pixelScale, //resolution
      maxPixels: 1e13
  });
}

processS2Data(rawImageStartDate, rawImageEndDate, selectBandList);
