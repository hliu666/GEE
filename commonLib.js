/***
 * 
 * get ablers china crs
 * */
function _getAlbersChinaCRS() {
  var wkt_text = 'PROJCS["unnamed",GEOGCS["WGS 84", DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]], AUTHORITY["EPSG","4326"]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["standard_parallel_1",25],PARAMETER["standard_parallel_2",47],PARAMETER["latitude_of_center",0],PARAMETER["longitude_of_center",105],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]';
  var albers_china = ee.Projection(wkt_text);
  return albers_china;
}

exports.getAlbersChinaCRS = function() {
  return _getAlbersChinaCRS();
};

/***
 * 
 * add gov shape outline
 * */
exports.addOutline = function(studyArea) {
  var empty = ee.Image().toByte();
  var outline = empty.paint({
    featureCollection:studyArea,
    color:0,
    width:3
  });
  Map.setOptions("SATELLITE");
  Map.centerObject(studyArea, 9);
  var shapeLayer = Map.addLayer(outline, {palette: "00ff00"}, "outline");
  return shapeLayer;
};

/***
 * 
 * calculate area
 * */
exports.calculateArea = function(img, key, studyArea, pixelScale) {
  var crs = _getAlbersChinaCRS();
  var countDict = img.reduceRegion({
      reducer: ee.Reducer.count(),
      geometry: studyArea,
      scale: pixelScale,
      maxPixels: 1e15,
      crs: crs
  });
  print(countDict);
  var pixelCount = ee.Dictionary(countDict).get(key);
  var area = ee.Number(pixelCount).multiply(pixelScale * pixelScale);
  print("Aera/m2 : " , area);
  print("Aera/mu: " , area.multiply(0.0015));
  print("Aera/acre : " , area.multiply(0.0001));
};

/***
 * 
 * export image to google drive
 * */
exports.exportImageToDrive = function(image, studyArea, pixelScale) {
    Export.image.toDrive({
        image: image,
        description: "enter description",
        // crs: "EPSG:4326",
        region: studyArea,
        scale: pixelScale,
        maxPixels: 1e13
    });
};

