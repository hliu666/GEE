var mcd12q2 = ee.ImageCollection("MODIS/006/MCD12Q2"),
    mod13a1 = ee.ImageCollection("MODIS/006/MOD13A1"),
    myd13a1 = ee.ImageCollection("MODIS/006/MYD13A1"),
    cfsv2 = ee.ImageCollection("NOAA/CFSV2/FOR6H"),
    countries = ee.FeatureCollection("ft:1tdSwUL7MVpOauSgRzqVTOwdfy17KDbw-1d9omPw"),
    geometry = table2; //Study Aera, Global, USA, China, etc


/**
1.Global SOS, EOS
MCD12Q2.005 Land Cover Dynamics Yearly Global 500m
SOS：Onset_Greenness_Increase1
EOS：Onset_Greenness_Minimum2

2. Avg NDVI、EVI、Meterological Variables in GSL
1）NDVI、EVI
Datasets：MOD13A1.006/MYD13A1.006

2）Temperature、Short-Wave、Precipitation (Temperature&Short-Wave: Avg, Precipitation:Accumulation)
Datasets：CFSV2: NCEP Climate Forecast System Version 2, 6-Hourly Products
Variables：Temperature_height_above_ground
Downward_Short-Wave_Radiation_Flux_surface_6_Hour_Average
Precipitation_rate_surface_6_Hour_Average

Calculate: GDD&EDD

I Linear Regression:
1.Detrend: set y as T、P、SW、GDD、EDD、NDVI、EVI、GSL; set X as year：y= a0+a1x

2.rawimage-simulated values = △T、△P、△SW、△GDD、△EDD、△NDVI、△EVI、△GSL

3.set △T、△P、△SW、△GDD、△EDD、△NDVI、△EVI、△GSL as x，set △NDVI、△EVI、△GSL as y seperately

**/
var region = countries.filter(ee.Filter.eq('Country', 'India'));
region = geometry;

var startDate = "2001-1-1";
var endDate = "2019-1-1";
//Filter time series are longer than 5 years
var yearNum = 5;
var selectKeys = [
  "Greenup_1",
  "Greenup_2",
  "Senescence_1",
  "Senescence_2",
  // "Onset_Greenness_Decrease1",
  // "Onset_Greenness_Decrease2",
  // "Onset_Greenness_Maximum1",
  // "Onset_Greenness_Maximum2"
];
var newSelectKeys = [
  "sday1",
  "sday2",
  "eday1",
  "eday2",
  "gsl",
  // "Onset_Greenness_Decrease1",
  // "Onset_Greenness_Decrease2",
  // "Onset_Greenness_Maximum1",
  // "Onset_Greenness_Maximum2"
];

//gdd and edd effect time
var t_base = 8;
var t_max = 30;
var t_edd = 34;
var rawLayer = null;
var pointLayer = null;
Map.centerObject(region, 4);
Map.style().set('cursor', 'crosshair');
// Map.addLayer(region, {color:"red"}, "region", false);

function maskImageByYearCount(imageCol, yearCount) {
  var meanImage = imageCol.mean();
  var sumImage = imageCol.sum();
  var mask1 = sumImage.select("Greenup_1")
                      .subtract(meanImage.select("Greenup_1").multiply(yearCount));
  mask1 = mask1.updateMask(mask1.gte(0));
  mask1 = mask1.unmask(0);
  var mask2 = sumImage.select("Greenup_2")
                      .subtract(meanImage.select("Greenup_2").multiply(yearCount));
  mask2 = mask2.updateMask(mask2.gte(0));
  mask2 = mask2.unmask(0);
  var mask3 = sumImage.select("Senescence_1")
                      .subtract(meanImage.select("Senescence_1").multiply(yearCount));
  mask3 = mask3.updateMask(mask3.gte(0));
  mask3 = mask3.unmask(0);
  var mask4 = sumImage.select("Senescence_2")
                      .subtract(meanImage.select("Senescence_2").multiply(yearCount));
  mask4 = mask4.updateMask(mask4.gte(0));
  mask4 = mask4.unmask(0);
  
  var mask = ee.Image.constant(0);
  mask = mask.where(mask1.or(mask2).and(mask3.or(mask4)), 1);
  
  imageCol = imageCol.map(function(image){
    return image.updateMask(mask);
  });
  return imageCol;
}

function processSeasonData() {
  var mcd12q2Col = mcd12q2.filterDate(startDate, endDate)
                          .map(function(image) {
                            var time_start = image.get("system:time_start");
                            var year = ee.Date(time_start).get("year");
                            var deltaStartDate = ee.Date("2000-1-1");
                            var deltaDay = ee.Date.fromYMD(year, 1, 1).difference(deltaStartDate, "day");
                            image = image.subtract(deltaDay).add(1);
                            image = image.set("year", year);
                            image = image.set("system:time_start", time_start);
                            return image.toInt();
                          })
                          .map(function(image){
                            return image.clip(region);
                          })
                          .select(selectKeys);
  mcd12q2Col = maskImageByYearCount(mcd12q2Col, yearNum);
  
  mcd12q2Col = mcd12q2Col.map(function(image){
    image = image.unmask(0);
    var ogi1 = image.select("Greenup_1");
    var ogi2 = image.select("Greenup_2");
    var ogm1 = image.select("Senescence_1");
    var ogm2 = image.select("Senescence_2");
    //I
    var mask1 = ogi1.gt(0).and(ogm1.gt(ogi1));
    var img1 = ogm1.subtract(ogi1);
    img1 = img1.updateMask(mask1)
               .unmask(0);
    var mask2 = ogi2.gt(0).and(ogm2.gt(ogi2));
    var img2 = ogm2.subtract(ogi2);
    img2 = img2.updateMask(mask2)
               .unmask(0);
    var delta1 = ogm1.subtract(ogi2);
    delta1 = delta1.updateMask(delta1.gt(0).and(ogi2.gt(0)))
                   .unmask(0);
    var gsl1 = img1.add(img2).subtract(delta1);
    //II
    var mask3 = ogm1.gt(0).and(ogi1.gt(ogm1));
    var img3 = ogm1.updateMask(mask3)
                   .unmask(0);
    var mask4 = ogm2.gt(0).and(ogi2.gt(ogm2));
    var yearImage = ee.Image.constant(365)
                      .updateMask(ogi2)
                      .unmask(0);
    var img4 = yearImage.subtract(ogi2);
    img4 = img4.updateMask(mask4)
               .unmask(0);
    var img5 = ogm2.subtract(ogi1);
    img5 = img5.updateMask(img5.gt(0).and(ogi1.gt(0)).and(mask3.or(mask4)))
               .unmask(0);
    var gsl2 = img3.add(img4).add(img5);
    //III
    var mask6 = ogi1.eq(0).and(ogi2.gt(0)).and(ogm1.gt(0)).and(ogm2.eq(0));
    var img6 = ogm1.add(yearImage.subtract(ogi2));
    var gsl3 = img6.updateMask(img6.gt(0).and(mask6))
                   .unmask(0);
    //VI
    var mask7 = ogi1.gt(0).and(ogi2.eq(0)).and(ogm1.eq(0)).and(ogm2.gt(0));
    var img7 = ogm2.subtract(ogi1);
    var gsl4 = img7.updateMask(img7.gt(0).and(mask7))
                   .unmask(0);
    
    var gsl = gsl1.add(gsl2).add(gsl3).add(gsl4);
    image = image.addBands(gsl.updateMask(gsl.gt(0)).rename("gsl"));
    image = image.addBands(ogi1.rename("sday1"));
    image = image.addBands(ogi2.rename("sday2"));
    image = image.addBands(ogm1.rename("eday1"));
    image = image.addBands(ogm2.rename("eday2"));
    return image;
  })
  .select(newSelectKeys);
  
  // addPanel(mcd12q2Col, "system:index");
  var meanGSL = mcd12q2Col.mean().toInt();
  var visParam = {min:1, max:365, palette:"#0f70a9,#8cc15e,#fee537,#d62631"};
  Map.addLayer(meanGSL.select("gsl").updateMask(meanGSL.select("gsl").gt(0)), visParam, "meanGSL");
  // clickMap(mcd12q2Col, newSelectKeys);
  return mcd12q2Col;
}

/**
 * 
 * @param {[type]} image [description]
 */
function setDeltaDay(image) {
  var doy = ee.Number.parse(ee.Date(image.get("system:time_start")).format("D"));
  image = image.addBands(ee.Image.constant(doy).rename("DayOfYear").toInt());
  return image;
}

/**
 * Mask images based on GSL
 * @param  {[type]} image [description]
 * @param  {[type]} sday1 [description]
 * @param  {[type]} eday1 [description]
 * @param  {[type]} sday2 [description]
 * @param  {[type]} eday2 [description]
 * @return {[type]}       [description]
 */
function maskImageBySeasonDay(image, sday1, eday1, sday2, eday2) {
  var deltaDay = image.select("DayOfYear");
  var mask11 = eday1.gte(sday1).and(sday1.gt(0))
                   .or(eday2.gte(sday2).and(sday2.gt(0)));
  var mask12 = deltaDay.gte(sday1).and(deltaDay.lte(eday1))
                       .or(deltaDay.gte(sday2).and(deltaDay.lte(eday2)));
  
  var mask21 = sday1.gte(eday1).and(eday1.gt(0))
                   .or(sday2.gte(eday2).and(eday2.gt(0)));
  var mask22 = deltaDay.lte(eday1)
                      .or(deltaDay.gte(sday2))
                      .or(deltaDay.gte(sday1).and(deltaDay.lte(eday2)));
                      
  var mask31 = sday1.eq(0)
                    .and(eday1.gt(0))
                    .and(sday2.gt(0))
                    .and(eday2.eq(0));
  var mask32 = deltaDay.lte(eday1).and(deltaDay.gte(sday2));
  
  
  var mask41 = sday1.gt(0)
                    .and(eday1.eq(0))
                    .and(sday2.eq(0))
                    .and(eday2.gt(0));
  var mask42 = deltaDay.lte(eday2).and(deltaDay.gte(sday1));
  
  
  var mask = mask11.and(mask12)
                   .or(mask21.and(mask22))
                   .or(mask31.and(mask32))
                   .or(mask41.and(mask42));
  image = image.updateMask(mask);
  return image;
}

function generateImageCol() {
  var mcd12q2Col = processSeasonData();
  
  var result = mcd12q2Col.map(function(mcdImg) {
    var year = ee.Number(mcdImg.get("year"));
    var sdate = ee.Date.fromYMD(year, 1, 1);
    var edate = ee.Date.fromYMD(year.add(1), 1, 1);
    
    var mod13a1Col = mod13a1.filterDate(sdate, edate);
    var myd13a1Col = myd13a1.filterDate(sdate, edate);
    var modisCol = mod13a1Col.merge(myd13a1Col);
    var cfsv2Col = cfsv2.filterDate(sdate, edate)
                        .map(setDeltaDay);
    
    var sday1 = mcdImg.select("sday1");
    var eday1 = mcdImg.select("eday1");
    var sday2 = mcdImg.select("sday2");
    var eday2 = mcdImg.select("eday2");
    ///////////////////////////////////////////////////////////////
    ///////////////year///////////////////////
    mcdImg = mcdImg.addBands(ee.Image.constant(year)
                               .toInt()
                               .rename("year"));
      
    ///////////////////////////////////////////////////////////////
    ///////////////gsl///////////////////////
    var gslImg = mcdImg.select("gsl");
    
    ///////////////////////////////////////////////////////////////
    ///////////////NDVI EVI///////////////////////
    var indexImg = modisCol.select(["NDVI", "EVI", "DayOfYear"])
                       .map(function(image) {
                         return maskImageBySeasonDay(image, sday1, eday1, sday2, eday2);
                       })
                       .mean()
                       .multiply(0.0001);
    mcdImg = mcdImg.addBands(indexImg.select("NDVI"));
    mcdImg = mcdImg.addBands(indexImg.select("EVI"));
    
    
    ///////////////////////////////////////////////////////////////
    var scol = cfsv2Col.select([
                          "Temperature_height_above_ground", 
                          "Downward_Short-Wave_Radiation_Flux_surface_6_Hour_Average",
                          "Precipitation_rate_surface_6_Hour_Average",
                          "DayOfYear"
                        ])
                        .map(function(image) {
                          return maskImageBySeasonDay(image, sday1, eday1, sday2, eday2);
                        });
                        
    var temCol = scol.select("Temperature_height_above_ground")
                      .map(function(image) {
                        image = image.subtract(273.15);
                        return image;
                      });
    ///////////////tem///////////////////////
    var temImg = temCol.mean();
    mcdImg = mcdImg.addBands(temImg.rename("tem"));
    
    ///////////////gdd///////////////////////
    var gddImg = temCol.map(function(image) {
                      var newImage = image.where(image.lt(t_base), 0);
                      newImage = newImage.where(newImage.gt(t_max), t_max - t_base);
                      image = image.addBands(newImage.rename("gdd"));
                      return image;
                    })
                    .sum()
                    .divide(4);
    mcdImg = mcdImg.addBands(gddImg.select("gdd"));
    
    ///////////////edd///////////////////////
    var eddImg = temCol.map(function(image) {
                      var newImage = image.where(image.lt(t_edd), 0);
                      image = image.addBands(newImage.rename("edd"));
                      return image;
                    })
                    .sum()
                    .divide(4);
    mcdImg = mcdImg.addBands(eddImg.select("edd"));
    
    ///////////////wave///////////////////////
    var waveImg = scol.select("Downward_Short-Wave_Radiation_Flux_surface_6_Hour_Average")
                      .mean();
    mcdImg = mcdImg.addBands(waveImg.rename("wave"));
    
    ///////////////pre///////////////////////
    var preImg = scol.select("Precipitation_rate_surface_6_Hour_Average")
                      .sum()
                      .multiply(6*60*60);
    mcdImg = mcdImg.addBands(preImg.rename("pre"));
    
    mcdImg = mcdImg.updateMask(gslImg.gt(0));
    mcdImg = mcdImg.clip(region);
    return mcdImg;
  });
  result = result.select(["gsl", "year", "NDVI", "EVI", "tem", "gdd", "edd", "wave", "pre"]);
  print("generateImageCol image collection is", result);
  return result;
}

function addVariables(image) {
  return image.addBands(ee.Image.constant(1).rename("constant"));
}

/**
 * Linear regression for imageCollection
 * @param  {[type]} imageCollection 
 * @param  {[type]} xValues         
 * @param  {[type]} yValue          
 * @param  {[type]} needDelta       
 * @return {[type]}                
 */
function linearRegressionFitted(imageCollection, xValues, yValue, needDelta) {
  //Add Constant
  var newImgCol = imageCollection.map(addVariables);
  xValues = ee.List(xValues);

  var trend = newImgCol.select(xValues.add(ee.String(yValue)))
                       .reduce(ee.Reducer.linearRegression(xValues.length(), 1));

  var coefficients = trend.select("coefficients")
                          .arrayProject([0])
                          .arrayFlatten([xValues]);
  if (!needDelta) {
    print(yValue + " coefficients is", coefficients);
  }

  var fittedImgCol = newImgCol.map(function(image) {

    image = xValues.iterate(function(xValue, img) {
      xValue = ee.String(xValue);
      img = ee.Image(img);
      var coef = coefficients.select(xValue);
      return img.addBands(coef.rename(xValue.cat("_c")));
    }, image);
    image = ee.Image(image);
    var fittedImg = image.select(xValues)
                        .multiply(coefficients)
                        .reduce("sum")
                        .rename("fitted");
    image = image.addBands(fittedImg);
    if (needDelta) {
      var deltaImg = image.select(yValue).subtract(fittedImg).rename("delta");
      image = image.addBands(deltaImg);
    }
    return image;
  });
  
  return fittedImgCol;
}

///////////////////////////////result/////////////////////////////

function generateDeltaImageCol() {
  var imgCol = generateImageCol();
  // TEM
  var temValues = ee.List(['constant', 'year']);
  var temValue = 'tem';
  var temFitCol = linearRegressionFitted(imgCol, temValues, temValue, true);
  
  // PRE
  var preValues = ee.List(['constant', 'year']);
  var preValue = 'pre';
  var preFitCol = linearRegressionFitted(imgCol, preValues, preValue, true);
  
  // WAVE
  var waveValues = ee.List(['constant', 'year']);
  var waveValue = 'wave';
  var waveFitCol = linearRegressionFitted(imgCol, waveValues, waveValue, true);
  
  // GDD
  var gddValues = ee.List(['constant', 'year']);
  var gddValue = 'gdd';
  var gddFitCol = linearRegressionFitted(imgCol, gddValues, gddValue, true);
  
  // EDD
  var eddValues = ee.List(['constant', 'year']);
  var eddValue = 'edd';
  var eddFitCol = linearRegressionFitted(imgCol, eddValues, eddValue, true);
  
  // NDVI
  var ndviValues = ee.List(['constant', 'year']);
  var ndviValue = 'NDVI';
  var ndviFitCol = linearRegressionFitted(imgCol, ndviValues, ndviValue, true);
  
  // EVI
  var eviValues = ee.List(['constant', 'year']);
  var eviValue = 'EVI';
  var eviFitCol = linearRegressionFitted(imgCol, eviValues, eviValue, true);
  
  // GSL
  var gslValues = ee.List(['constant', 'year']);
  var gslValue = 'gsl';
  var gslFitCol = linearRegressionFitted(imgCol, gslValues, gslValue, true);
  
  
  var deltaImgCol = imgCol.select(["year", "gsl"])
                          .combine(temFitCol.select(["delta"], ["deltaTEM"]))
                          .combine(preFitCol.select(["delta"], ["deltaPRE"]))
                          .combine(waveFitCol.select(["delta"], ["deltaWAVE"]))
                          .combine(gddFitCol.select(["delta"], ["deltaGDD"]))
                          .combine(eddFitCol.select(["delta"], ["deltaEDD"]))
                          .combine(ndviFitCol.select(["delta"], ["deltaNDVI"]))
                          .combine(eviFitCol.select(["delta"], ["deltaEVI"]))
                          .combine(gslFitCol.select(["delta"], ["deltaGSL"]));
  print("deltaImgCol", deltaImgCol);
  return deltaImgCol;
}

/**
 * 添加图层到地图上
 * @param {[type]} imgCol [description]
 * @param {[type]} xKeys  [description]
 * @param {[type]} yKey   [description]
 * @param {[type]} flag   [description] 是否显示结果图层
 */
function addLayerToMap(imgCol, xKeys, yKey, flag) {
  var newXKeys = [];
  for (var j=0; j<xKeys.length; j++) {
    newXKeys.push(xKeys[j]+"_c");
  }
  var coefImg = ee.Image(imgCol.first()).select(newXKeys);
  //export image
  exportResult(coefImg, ee.FeatureCollection(region).geometry().bounds(), yKey); 
  
  //show layer
  if (flag) {
    var visParam = {
      min:-0.1,
      max:0.1,
      palette: 'FFFFFF, CE7E45, DF923D, F1B555, FCD163, 99B718, 74A901, 66A000, 529400,'+ 
               '3E8601, 207401, 056201, 004C00, 023B01, 012E01, 011D01, 011301'
    };
    for (var i = 0; i<newXKeys.length; i++) {
      Map.addLayer(coefImg.select(newXKeys[i]), visParam, yKey + "_" + newXKeys[i]);
    }
  }
}

/**
 * @param  {[type]} image  [description]
 * @param  {[type]} region [description]
 * @param  {[type]} key    [description]
 * @return {[type]}        [description]
 */
function exportResult(image, region, key) {
    Export.image.toDrive({
        image: image,
        description: key,
        crs: "EPSG:4326",
        region: region,
        scale: 500,
        maxPixels: 1e13
    });
}

function cacluateResult() {
  var deltaImgCol = generateDeltaImageCol();
  var deltaNDVIValue = 'deltaNDVI';
  var deltaGSLValue = 'deltaGSL';
  var deltaEVIValue = 'deltaEVI';
  
  //Reg1: Y:△NDVI/△EVI/△GSL， X:△T、△P、△WAVE
  var deltaValues1 = ['constant', 'deltaTEM', 'deltaPRE', 'deltaWAVE'];
 
  ////△NDVI
  //var deltaNDVIFitCol1 = linearRegressionFitted(deltaImgCol, deltaValues1, deltaNDVIValue, false);
  //print("deltaNDVIFitCol1", deltaNDVIFitCol1);
  //addLayerToMap(deltaNDVIFitCol1, deltaValues1, deltaNDVIValue+"_1", false);
  

  // //△EVI
  // var deltaEVIFitCol1 = linearRegressionFitted(deltaImgCol, deltaValues1, deltaEVIValue, false);
  // print("deltaEVIFitCol1", deltaEVIFitCol1);
  // addLayerToMap(deltaEVIFitCol1, deltaValues1, deltaEVIValue+"_1", false);
  
  // //△GSL
  // var deltaGSLFitCol1 = linearRegressionFitted(deltaImgCol, deltaValues1, deltaGSLValue, false);
  // print("deltaGSLFitCol1", deltaGSLFitCol1);
  // addLayerToMap(deltaGSLFitCol1, deltaValues1, deltaGSLValue+"_1", false);
  
  //Reg2: Y:△NDVI/△EVI/△GSL，X:△GDD、△EDD、△P、△WAVE
  var deltaValues2 = ['constant', 'deltaGDD', 'deltaEDD', 'deltaPRE', 'deltaWAVE'];
  
  // //△NDVI
  // var deltaNDVIFitCol2 = linearRegressionFitted(deltaImgCol, deltaValues2, deltaNDVIValue, false);
  // print("deltaNDVIFitCol2", deltaNDVIFitCol2);
  // addLayerToMap(deltaNDVIFitCol2, deltaValues2, deltaNDVIValue+"_2", false);
  
  // //△EVI
  // var deltaEVIFitCol2 = linearRegressionFitted(deltaImgCol, deltaValues2, deltaEVIValue, false);
  // print("deltaEVIFitCol2", deltaEVIFitCol2);
  // addLayerToMap(deltaEVIFitCol2, deltaValues2, deltaEVIValue+"_2", false);
  
   // △GSL
   var deltaGSLFitCol2 = linearRegressionFitted(deltaImgCol, deltaValues2, deltaGSLValue, false);
   print("deltaGSLFitCol2", deltaGSLFitCol2);
   addLayerToMap(deltaGSLFitCol2, deltaValues2, deltaGSLValue+"_2", true);
  
  // Map.onClick(function(coords) {
  //   print("click map point is: " + coords.lon + " " + coords.lat);
  //   var point = ee.Geometry.Point(coords.lon, coords.lat);
  //   showPointFittedData(deltaGSLFitCol1, point, deltaValues1, ["fitted", deltaGSLValue]);
  // });
}

//main
cacluateResult();

///////////////////////////////////////////// TEST /////////////////////////////////
function addPanel(scol, keyName) {
  var id_list = scol.reduceColumns(ee.Reducer.toList(), [keyName])
                    .get('list');
  print("id list is", id_list);
  id_list.evaluate(function(ids) {
    var total = ids.length;
    var showTitle = ui.Label("", {fontWeight: 'bold'});
    var curIndex = 0;
    var bPlus = ui.Button("+", function() {
      curIndex += 1;
      if (curIndex >= total) {
        curIndex = 0;
      }
      showTitle.setValue(ids[curIndex]);
      showSelectRawImage(scol, keyName, ids[curIndex]);
    });
    var bReduce = ui.Button("-", function() {
      curIndex -= 1;
      if (curIndex < 0) {
        curIndex = total - 1;
      }
      showTitle.setValue(ids[curIndex]);
      showSelectRawImage(scol, keyName, ids[curIndex]);
    });
    showTitle.setValue(ids[curIndex]);
    showSelectRawImage(scol, keyName, ids[curIndex]);
    var main = ui.Panel({
      widgets: [
        ui.Label('click "+" or "-" to move time window', {fontWeight: 'bold'}),
        bPlus, bReduce,
        ui.Label("select date: ", {fontWeight: 'bold'}),
        showTitle
      ],
      style: {width: '210px', padding: '8px'}
    });
    ui.root.insert(0, main);
  });
}

function showSelectRawImage(scol, keyName, keyValue) {
  if (rawLayer !== null) {
    Map.remove(rawLayer);
  }
  print("show raw image id is: " + keyValue);
  var simg = scol.filter(ee.Filter.eq(keyName, keyValue)).first();
  rawLayer = Map.addLayer(
    simg.updateMask(simg.select("gsl")),
    {bands:["gsl"], min:1, max:365, palette:"blue,green,red"},
    "gsl"
  );
}

function clickMap(imgCol, showKeys) {
  Map.onClick(function(coords) {
    if (pointLayer !== null) {
      Map.remove(pointLayer);
    }
    print("click map point is: " + coords.lon + " " + coords.lat);
    var point = ee.Geometry.Point(coords.lon, coords.lat);
    pointLayer = Map.addLayer(point, {color:"white"}, "point");
    
    var data_list = imgCol.getRegion(point, 500);
    data_list.evaluate(function(result){
      var band_value_list = [];
      for(var i=1; i<result.length; i++) {
        var temp = {};
        for (var j=0; j<showKeys.length; j++) {
          temp[showKeys[j]] = result[i][j+4];
        }
        band_value_list.push(temp);
      }
      print("season data list", band_value_list);
    });
  });
}

