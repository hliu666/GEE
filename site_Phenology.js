var mod13a2 = ee.ImageCollection("MODIS/006/MOD13A2"),
    myd13a2 = ee.ImageCollection("MODIS/006/MYD13A2"),
    province = ee.FeatureCollection("***/HeiLongJiang"),
    agri_stations = ee.FeatureCollection("***/agri_sites_utf8");

 
var startYear = 2007;
var endYear = 2015;
//xxxx 230000
//var roi = province.filter(ee.Filter.eq("ad2004", 230000));
var roi = province.filter(ee.Filter.eq("ADMINCODE", "232722"));
var pixelScale = 1000;

var modCol = mod13a2.filterDate(ee.Date.fromYMD(startYear, 1, 1), ee.Date.fromYMD(endYear,1,1))
                    .map(function(img) {
                      var time_start = img.get("system:time_start");
                      img = img.updateMask(img.select("SummaryQA").neq(3));
                      img = img.unitScale(0, 10000);
                      img = img.set("system:time_start", time_start);
                      return img;
                    }).select(["NDVI", "EVI"]);
                    
var mydCol = myd13a2.filterDate(ee.Date.fromYMD(startYear, 1, 1), ee.Date.fromYMD(endYear,1,1))
                    .map(function(img) {
                      var time_start = img.get("system:time_start");
                      img = img.updateMask(img.select("SummaryQA").neq(3));
                      img = img.unitScale(0, 10000);
                      img = img.set("system:time_start", time_start);
                      return img;
                    }).select(["NDVI", "EVI"]);
                    
var modisCol = modCol.merge(mydCol).sort("system:time_start");
print(modisCol.limit(10));

var visParam = {
  min: -0.2, 
  max: 0.8,
  palette: 'FFFFFF, CE7E45, DF923D, F1B555, FCD163, 99B718, 74A901, 66A000, 529400,' +
    '3E8601, 207401, 056201, 004C00, 023B01, 012E01, 011D01, 011301'
};
Map.addLayer(ee.Image(mydCol.first()).clip(roi).select("NDVI"), visParam, "NDVI");

///////////////////////////////////////////////////// 
var SG_Models = {
  /***
  * 
  * s-g filter
  * */
  sgFilter : function(y, window_size, order) {
    var half_window = (window_size - 1)/2;
    var deriv = 0;
    var order_range = ee.List.sequence(0,order);
    var k_range = ee.List.sequence(-half_window, half_window);
    
    //b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    var b = ee.Array(k_range.map(function (k) { return order_range.map(function(o) { return ee.Number(k).pow(o)})}));
    // m = np.linalg.pinv(b).A[deriv] 
    var mPI = b.matrixPseudoInverse();
    var impulse_response = (mPI.slice({axis: 0, start: deriv, end: deriv+1})).project([1]);
    //firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    var y0 = y.get(0);
    var firstvals = y.slice(1, half_window+1).reverse().map(
      function(e) { return ee.Number(e).subtract(y0).abs().multiply(-1).add(y0) }
    );
    
    //lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    var yend = y.get(-1);
    var lastvals = y.slice(-half_window-1,-1).reverse().map(
      function(e) { return ee.Number(e).subtract(yend).abs().add(yend) }
    );
    
    // y = np.concatenate((firstvals, y, lastvals))
    var y_ext = firstvals.cat(y).cat(lastvals);
    
    // np.convolve( m, y, mode='valid')
    var runLength = ee.List.sequence(0, y_ext.length().subtract(window_size));
    
    var smooth = runLength.map(function(i) {
      return ee.Array(y_ext.slice(ee.Number(i), ee.Number(i).add(window_size))).multiply(impulse_response).reduce("sum", [0]).get([0]);
    });
    return smooth;
  },
  
  /***
  * 
  * fill nodata
  * */
  fill_nodata: function(index, data_list, step) {
    var result = null;
    var i1 = (index - step >= 0) ? (index - step) : 0;
    var i2 = (index + step <= data_list.length - 1) ? (index + step) : data_list.length - 1;
    var max1 = 0;
    var max2 = 0;
    var i = 0;
    for (i=i1; i<index; i++) {
      if(data_list[i] !== null) {
        if (max1 < data_list[i]) {
          max1 = data_list[i];
        }
      }
    }
    
    i = 0;
    for (i=index; i<i2; i++) {
      if(data_list[i] !== null) {
        if (max2 < data_list[i]) {
          max2 = data_list[i];
        }
      }
    }
    result = (max1 + max2)/2;
    return result;
  },
  
  addMultiSmoothImageDoySeriesByPoint: function(roi, dataset, bandNames, pixelScale) {
    var raw_data_list = dataset.select(bandNames).getRegion(roi, pixelScale);
    var bandCount = bandNames.length;
    raw_data_list.evaluate(function(result) {
      var bands_value_list = [];
      for(var j=0; j<bandCount; j++) {
        bands_value_list.push([]);
      }
      
      var time_list = [];
      for (var i =0; i<result.length; i++) {
        var data = result[i];
        if (data[0] !== "id") {
          var time_str = ee.Date(data[3]).format("YYYY-MM-dd");
          time_list.push(time_str);
          for(var k=0; k<bandCount; k++) {
            bands_value_list[k].push(data[4+k]);
          }
        }
      }
      
      var full_data_list = [];
      for(j=0; j<bandCount; j++) {
        full_data_list.push([]);
      }
      for (i = 0; i < bands_value_list.length; i++) {
        var d_list = bands_value_list[i];
        for (j=0; j<d_list.length; j++) {
          var d = d_list[j];
          d = d === null ? SG_Models.fill_nodata(j, d_list, 3): d;
          full_data_list[i].push(d);
        }
      }
      
      var sg_filter_list = [];
      var window_size = 7;
      var order = 2;
      for (i=0; i<full_data_list.length; i++) {
        var src_list = ee.List(full_data_list[i]);
        var smooth = SG_Models.sgFilter(src_list, window_size, order);
        sg_filter_list.push(smooth);
      }
      var yValues = ee.Array.cat(sg_filter_list, 1);
      print("sg filter list :", yValues);
      var chart = ui.Chart.array.values(yValues, 0, ee.List(time_list))
                    .setSeriesNames(bandNames)
                    .setOptions({
                      title: 'SG-Model Window: ' + window_size + ' Order: ' + order, 
                      hAxis: {title: 'time'}, vAxis: {title: "Smooth Value"},
                      legend: null,
                      lineWidth:1,
                      pointSize:2
                    });
      print(chart);
    });
  },

  addSingleSmoothImageDoySeriesByPoint: function(roi, dataset, bandName, pixelScale) {
    //smooth
    var raw_data_list = dataset.select(bandName).getRegion(roi, pixelScale);
    // print("raw_data_list is : ", raw_data_list);
    
    raw_data_list.evaluate(function(result) {
      var band_value_list = [];
      var times_list = [];
      for(var j=1; j<result.length; j++) {
        var data = result[j];
        var time_str = ee.Date(data[3]).format("YYYY-MM-dd");
        times_list.push(time_str);
        band_value_list.push(data[4]);
      }
      
      var full_data_list = [];
      for (var i = 0; i < band_value_list.length; i++) {
        var d = band_value_list[i];
        if (d === null) {
          d = SG_Models.fill_nodata(i, band_value_list, 3);
        }
        full_data_list.push(d);
      }
      var window_size = 7;
      var order = 2;
      var smooth_list = SG_Models.sgFilter(ee.List(full_data_list), window_size, order);
      print("smooth_list is: ", smooth_list);
      var chart = ui.Chart.array.values(smooth_list, 0, ee.List(times_list))
                    .setSeriesNames(['Smoothed'])
                    .setOptions({
                      title: 'SG-Model Window: ' + window_size + ' Order: ' + order, 
                      hAxis: {title: 'time'}, vAxis: {title: bandName},
                      legend: null,
                      series: { 
                        0: { lineWidth: 1, pointSize: 2 }
                      }
                    });
      print(chart);
    });
  }
};

function addSingleSmoothImageDoySeriesByPoint(roi, dataset, bandName, pixelScale) {
  //smooth
  var raw_data_list = dataset.select(bandName).getRegion(roi, pixelScale);
  raw_data_list = raw_data_list.slice(1);
  raw_data_list.evaluate(function(result) {
    var band_value_list = [];
    var date_list = [];
    for(var j=0; j<result.length; j++) {
      var data = result[j];
      var date_str = ee.Date(data[3]).format("YYYY-MM-dd");
      date_list.push(date_str);
      band_value_list.push(data[4]);
    }
    date_list = ee.List(date_list);
    
    var full_data_list = [];
    for (var i = 0; i < band_value_list.length; i++) {
      var d = band_value_list[i];
      if (d === null) {
        d = SG_Models.fill_nodata(i, band_value_list, 3);
      }
      full_data_list.push(d);
    }
    var window_size = 7;
    var order = 2;
    var smooth_list = SG_Models.sgFilter(ee.List(full_data_list), window_size, order);
    // print("smooth_list is: ", smooth_list);
    
    var mean_value = getSeasonBaseValue(smooth_list);
    var season_data_list = getSeasonDataList(date_list, smooth_list, mean_value);
    var max_value_list = getSeasonMaxValueList(season_data_list);
    var max_value = getSeasonMaxValue(max_value_list);
    var season_param_list = calculateSeasonParams(season_data_list, max_value_list, mean_value);
    var season_date_list = calculateSeasonDateList(date_list, smooth_list, season_param_list);
    
    var y_list = ee.List(smooth_list.iterate(function(data, list){
      var temp = ee.List([
        ee.Number(data), max_value, mean_value, season_date_list.get(ee.List(list).size())
      ]);
      return ee.List(list).add(temp);
    }, ee.List([])));
    var x_list = date_list;
    var chart = ui.Chart.array.values(y_list, 0, x_list)
                  .setSeriesNames(['Smoothed ' + bandName , 'MaxValue','BaseValue', 'Season Date'])
                  .setOptions({
                    title: 'SG-Model Window: ' + window_size + ' Order: ' + order, 
                    hAxis: {title: 'time'}, vAxis: {title: bandName},
                    legend: null,
                    series: { 
                      0: { lineWidth: 1, pointSize: 1.5},
                      1: { lineWidth: 1, pointSize: 0, color: "00FFFF"},
                      2: { lineWidth: 1, pointSize: 0, color: "00FF00"},
                      3: { lineWidth: 0, pointSize: 3, color: "FF0000"}
                    }
                  });
    print(chart);
  });
}

/***
 * calculate detrended
 * */
function getDetrendedDataList(date_list, smooth_list) {
  var a_list = smooth_list.slice(0, ee.Number(smooth_list.size()).subtract(1));
  var b_list = smooth_list.slice(1);
  var c_list = a_list.zip(b_list);
  var detrended_list = c_list.iterate(function(data, list){
    var data_list = ee.List(data);
    var detrend = ee.Number(data_list.get(1)).subtract(ee.Number(data_list.get(0)));
    return ee.List(list).add(detrend);
  }, ee.List([]));
  detrended_list = ee.List(detrended_list);
  detrended_list = detrended_list.add(0);
  print("detrended_list", detrended_list);
  
  var chart = ui.Chart.array.values(detrended_list, 0, date_list)
                .setSeriesNames(['Detrended Value'])
                .setOptions({
                  title: 'Detrended Value list', 
                  hAxis: {title: 'time'}, vAxis: {title: "value"},
                  legend: null,
                  lineWidth: 1, 
                  pointSize: 2
                });
  print(chart);
  return detrended_list;
}

/**
 * This amplitude is calculated as the diﬀerence between the robust mean maximum and the robust mean base level
 * (the means of values when excluding the 10 % lowest and highest values) 
 **/
function getSeasonBaseValue(smooth_list) {
  var percentile = smooth_list.reduce(ee.Reducer.percentile([10, 90]));
  percentile = ee.Dictionary(percentile);
  // print("getSeasonBaseValue percentile", percentile);
  var p10 = percentile.get("p10");
  var p90 = percentile.get("p90");
  var select_data_list = smooth_list.map(function(data) {
    data = ee.Number(data);
    var new_data = ee.Algorithms.If(data.gt(p90).or(data.lt(p10)), null, data);
    return new_data;
  });
  var mean_value = select_data_list.reduce(ee.Reducer.mean());
  print("mean_value is: ", mean_value);
  return mean_value;
}

/***
 * calcluate season data list
 * */
function getSeasonDataList(date_list, smooth_list, mean_value) {
  var smooth_date_value_list = date_list.zip(smooth_list);
  // print("getSeasonDataList smooth_date_value_list", smooth_date_value_list);
  
  //get all great base value data list
  var temp_season_data_list = smooth_date_value_list.iterate(function(data, list){
    data = ee.List(data);
    list = ee.List(list);
    var _date = ee.Date(data.get(0));
    var _value = ee.Number(data.get(1));
    list = ee.Algorithms.If(_value.gte(mean_value), list.add(data), list);
    return list;
  }, ee.List([]));
  temp_season_data_list = ee.List(temp_season_data_list);
  // print("getSeasonDataList temp_season_data_list", temp_season_data_list);
  
  var func1 = function(data, list) {
    var temp = list.add(ee.List([data]));
    return temp;
  };
  var func2 = function(data, list) {
    var _date = ee.Date(data.get(0));
    var pre_list = ee.List(list.get(-1));
    var pre_data = ee.List(pre_list.get(-1));
    var _pre_date = ee.Date(pre_data.get(0));
    var condition = ee.Number(_date.difference(_pre_date, "day")).gt(ee.Number(8));
    var temp = ee.Algorithms.If(condition, list.add(ee.List([data])), list.set(-1, pre_list.add(data)));
    return temp;
  };
  //get season data list
  var season_data_list = temp_season_data_list.iterate(function(data, list){
    data = ee.List(data);
    list = ee.List(list);
    var new_list = ee.Algorithms.If(list.size(), func2(data, list), func1(data, list));
    return new_list;
  }, ee.List([]));
  season_data_list = ee.List(season_data_list);
  print("getSeasonDataList season_data_list", season_data_list);
  return season_data_list;
}

/***
 * calcluate season max value list
 * */
function getSeasonMaxValueList(season_data_list) {
  var max_value_list = season_data_list.iterate(function(data_list, list){
    data_list = ee.List(data_list);
    var new_data_list = data_list.iterate(function(data_list, list){
      data_list = ee.List(data_list);
      list = ee.List(list).add(data_list.get(1));
      return list;
    }, ee.List([]));
    var max_value = ee.List(new_data_list).reduce(ee.Reducer.max());
    list = ee.List(list).add(max_value);
    return list;
  }, ee.List([]));
  max_value_list = ee.List(max_value_list);
  print("getSeasonMaxValueList max_value_list", max_value_list);
  return max_value_list;
}

/***
 * calculate season max value
 * */
function getSeasonMaxValue(max_value_list) {
  max_value_list = ee.List(max_value_list);
  var max_mean_value = max_value_list.reduce(ee.Reducer.mean());
  print("getSeasonMaxValue max_mean_value", max_mean_value);
  return max_mean_value;
}

/**
 * calcluate season amplitude
 * */
function calcluateAmplitude(max_value, mean_value) {
  var amplitude = ee.Number(max_value).subtract(ee.Number(mean_value));
  print("calcluateAmplitude season amplitude is: ", amplitude);
  return amplitude;
}

/***
 * calcluate differ season params
 * */
function calculateSeasonParams(season_data_list, max_value_list, mean_value) {
  var season_num = season_data_list.size();
  print("calculateSeasonParams season_num is", season_num);
  
  var season_param_list = season_data_list.iterate(function(data_list, list){
    var max_value = ee.Number(max_value_list.get(ee.List(list).size()));
    data_list = ee.List(data_list);
    var data_value_list = data_list.iterate(function(data_list, list){
      return ee.List(list).add(ee.List(data_list).get(1));
    }, ee.List([]));
    var max_index = ee.List(data_value_list).indexOf(max_value);
    var left_data = data_list.slice(0, ee.Number(max_index).add(1));
    var left_value_list = left_data.iterate(function(data, list){
      var _value = ee.Number(ee.List(data).get(1));
      return ee.List(list).add(_value);
    }, ee.List([]));
    var left_per = ee.List(left_value_list).reduce(ee.Reducer.percentile([10, 90]));
    var left_p10 = ee.Dictionary(left_per).get("p10");
    var left_p10_data = ee.List(left_data).iterate(function(data_list, list){
      data_list = ee.List(data_list);
      return ee.Algorithms.If(
        list, 
        list, 
        ee.Algorithms.If(
          ee.Number(data_list.get(1)).gte(ee.Number(left_p10)), 
          data_list, 
          ee.List([])
        )
      );
    }, ee.List([]));
    left_p10_data = ee.List(left_p10_data);
    
    var right_data = data_list.slice(ee.Number(max_index));
    var right_value_list = right_data.iterate(function(data, list){
      var _value = ee.Number(ee.List(data).get(1));
      return ee.List(list).add(_value);
    }, ee.List([]));
    var right_per = ee.List(right_value_list).reduce(ee.Reducer.percentile([10, 90]));
    var right_p10 = ee.Dictionary(right_per).get("p10");
    var right_p10_data = ee.List(right_data).iterate(function(data_list, list){
      data_list = ee.List(data_list);
      return ee.Algorithms.If(
        list, 
        list, 
        ee.Algorithms.If(
          ee.Number(data_list.get(1)).lte(ee.Number(right_p10)), 
          data_list, 
          ee.List([])
        )
      );
    }, ee.List([]));
    right_p10_data = ee.List(right_p10_data);
    
    // list = ee.List(list).add([left_value_list, left_p10]);
    var beginOfSeason = left_p10_data.get(0);
    var beginOfSeasonValue = ee.Number(left_p10_data.get(1));
    var endOfSeason = right_p10_data.get(0);
    var endOfSeasonValue = ee.Number(right_p10_data.get(1));
    var legthOfSeason = ee.Date(endOfSeason).difference(ee.Date(beginOfSeason), "day");
    var amplitude = max_value.subtract(mean_value);
    list = ee.List(list).add(ee.Dictionary({
      beginOfSeason: beginOfSeason,
      beginOfSeasonValue: beginOfSeasonValue,
      endOfSeason: endOfSeason,
      endOfSeasonValue: endOfSeasonValue,
      legthOfSeason: legthOfSeason,
      maxValue: max_value,
      baseValue: mean_value, 
      amplitude: amplitude
    }));
    return list;
  }, ee.List([]));
  season_param_list = ee.List(season_param_list);
  print("calculateSeasonParams season_param_list", season_param_list);
  return season_param_list;
}

function calculateSeasonDateList(date_list, smooth_list, season_param_list) {
  var select_date_list = season_param_list.iterate(function(data, list){
    data = ee.Dictionary(data);
    list = ee.List(list).add(data.get("beginOfSeason"));
    list = list.add(data.get("endOfSeason"));
    return list;
  }, ee.List([]));
  select_date_list = ee.List(select_date_list);
  print("calculateSeasonDateList select_date_list", select_date_list);
  
  var smooth_date_value_list = date_list.zip(smooth_list);
  var season_date_list = smooth_date_value_list.map(function(data){
    data = ee.List(data);
    var _date = data.get(0);
    var temp = ee.Algorithms.If(select_date_list.contains(_date), data.get(1), null);
    return temp;
  });
  season_date_list = ee.List(season_date_list);
  // print("calculateSeasonDateList season_date_list", season_date_list);
  return season_date_list;
}


function clickMap(point) {
  //normal
  var chart = ui.Chart.image.series({
      imageCollection: modisCol.select(["NDVI"]), 
      region:point, 
      reducer: ee.Reducer.mean(), 
      scale:pixelScale
  }).setOptions({
    title: 'Raw Data',
    vAxis: {title: 'Normal Value'},
    lineWidth: 1,
    pointSize: 2
  });
  print(chart);
  
  //smooth
  // SG_Models.addMultiSmoothImageDoySeriesByPoint(point, modisCol, ["EVI"], pixelScale);
  addSingleSmoothImageDoySeriesByPoint(point, modisCol, "NDVI", pixelScale);
}

function addBounds(area, zoom) {
  var empty = ee.Image().toByte();
  var bounds = empty.paint({
    featureCollection:area,
    color:0,
    width:1.5
  });
  Map.setOptions("SATELLITE");
  Map.centerObject(area, zoom);
  Map.style().set('cursor', 'crosshair');
  var shapeLayer = Map.addLayer(bounds, {palette: "0000ff"}, "bounds");
  return shapeLayer;
}

addBounds(roi, 7);
//click map
Map.onClick(function(coords) {
  print("click map point is: " + coords.lon + " " + coords.lat);
  var point = ee.Geometry.Point(coords.lon, coords.lat);
  clickMap(point);
});

//////////////////////////Site Data////////////////////////////////////////////
var show_stations = agri_stations.filterBounds(roi);
Map.addLayer(show_stations, {color:'red'}, "stations");
var data_list = show_stations.reduceColumns(ee.Reducer.toList(4), ["scode","sname", "lon", "lat"]).get("list");
// select check data
var codes = [50564,50756,50851,50964,50880,50892,54080];
var select_stations = agri_stations.filter(ee.Filter.inList("sname", ee.List(codes)));
Map.addLayer(select_stations, {color:"white"}, "select_stations");
