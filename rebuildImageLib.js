var s2col = ee.ImageCollection("COPERNICUS/S2");

 function mosaicByTime(rawImageCollection) {
  var TIME_FIELD = 'system:time_start';
  var distinct = rawImageCollection.distinct([TIME_FIELD]);
  
  var filter = ee.Filter.equals({ leftField: TIME_FIELD, rightField: TIME_FIELD });
  var join = ee.Join.saveAll('matches');
  var results = join.apply(distinct, rawImageCollection, filter);

  // mosaic
  results = results.map(function(i) {
    var mosaic = ee.ImageCollection.fromImages(i.get('matches')).sort('system:index').mosaic();
    return mosaic.copyProperties(i).set(TIME_FIELD, i.get(TIME_FIELD));
  });
  return ee.ImageCollection(results);
}

//S2
var S2 = {
  //remove cloud
  rmCloud: function(image) {
    var quality = image.select("QA60").unmask();
    return image.updateMask(quality.eq(0));
  },
  
  //NDBI: (B11 - B08)/(B11 + B08)
  NDBI: function(image) {
    return image.addBands(image.normalizedDifference(["B11", "B8"])
          .rename("NDBI"));
  },
  
  //NDWI: (B03 - B08)/(B03 + B08)
  NDWI: function(image) {
      return image.addBands(image.normalizedDifference(["B3", "B8"])
          .rename("NDWI"));
  },
  
  //NDVI: (B08 - B04)/(B08 + B04)
  NDVI: function(img) {
    var ndvi = img.normalizedDifference(["B8","B4"]);
    return img.addBands(ndvi.rename("NDVI"));
  },
  
  //EVI: 2.5*(B08 - B04) / (B08 + 6*B04 - 7.5*B02 + 1)
  EVI: function(img) {
    var nir = img.select("B8");
    var red = img.select("B4");
    var blue = img.select("B2");
    var evi = img.expression(
      "2.5 * (B8 - B4) / (B8 + 6*B4 - 7.5*B2 + 1)",
      {
        "B8": nir,
        "B4": red,
        "B2": blue
      }
    );
    return img.addBands(evi.rename("EVI"));
  },
  
  //EVI2: 2.5 * (B08 - B04) / (B08 + 2.4 * B04 + 1)
  EVI2: function(img) {
    var nir = img.select("B8");
    var red = img.select("B4");
    var evi2 = img.expression(
      "2.5 * (B8 - B4) / (B8 + 2.4 * B4 + 1)",
      {
        "B8": nir,
        "B4": red
      }
    );
    return img.addBands(evi2.rename("EVI2"));
  },
  
  //LSWI: (B08 - B11)/(B08 + B11)
  LSWI: function(img) {
    var lswi = img.normalizedDifference(["B8","B11"]);
    return img.addBands(lswi.rename("LSWI"));
  },
  
  //RCLE (LSWI_max - LSWI_min)/(EVI2_heading - EVI2_tilering)
  RCLE: function(img) {
    var LSWI_max = img.select("LSWI_max");
    var LSWI_min = img.select("LSWI_min");
    var EVI2_heading = img.select("EVI2_heading");
    var EVI2_tilering = img.select("EVI2_tilering");
  
    var rcle = img.expression(
      "(LSWI_max - LSWI_min)/(EVI2_heading - EVI2_tilering)",
      {
        "LSWI_max": LSWI_max,
        "LSWI_min": LSWI_min,
        "EVI2_heading": EVI2_heading,
        "EVI2_tilering": EVI2_tilering
      }
    );
    return img.addBands(rcle.rename("RCLE"));
  },
  
  filterS2Data : function(startDate, endDate, bounds) {
    var dataset = s2col.filterDate(startDate, endDate)
                      .filterBounds(bounds)
                      .select(["B2","B3","B4","B8","B11","QA60"])
                      .map(S2.rmCloud)
                      .map(S2.NDVI)
                      .map(S2.NDWI)
                      .map(S2.EVI)
                      .map(S2.EVI2)
                      .map(S2.LSWI)
                      .map(S2.NDBI)
                      .map(function(img){
                        return img.clip(bounds);
                      })
                      .sort("system:time_start");
    dataset = mosaicByTime(dataset);
    return dataset;
  },
  
  filterS2RawData : function(startDate, endDate, bounds) {
    var dataset = s2col.filterDate(startDate, endDate)
                      .filterBounds(bounds)
                      .select(["B2","B3","B4","B8","B11","QA60"])
                      .map(S2.rmCloud)
                      .map(function(img){
                        return img.clip(bounds);
                      })
                      .sort("system:time_start");
    dataset = mosaicByTime(dataset);
    return dataset;
  },
  
};


var SG_Models = {

  sgFilter : function(y, window_size, order) {
    var half_window = (window_size - 1)/2;
    var deriv = 0;
    var order_range = ee.List.sequence(0,order);
    var k_range = ee.List.sequence(-half_window, half_window);
    
    var b = ee.Array(k_range.map(function (k) { return order_range.map(function(o) { return ee.Number(k).pow(o)})}));
    var mPI = b.matrixPseudoInverse();
    var impulse_response = (mPI.slice({axis: 0, start: deriv, end: deriv+1})).project([1]);
    var y0 = y.get(0);
    var firstvals = y.slice(1, half_window+1).reverse().map(
      function(e) { return ee.Number(e).subtract(y0).abs().multiply(-1).add(y0) }
    );
    
    var yend = y.get(-1);
    var lastvals = y.slice(-half_window-1,-1).reverse().map(
      function(e) { return ee.Number(e).subtract(yend).abs().add(yend) }
    );
    
    var y_ext = firstvals.cat(y).cat(lastvals);
    var runLength = ee.List.sequence(0, y_ext.length().subtract(window_size));
    
    var smooth = runLength.map(function(i) {
      return ee.Array(y_ext.slice(ee.Number(i), ee.Number(i).add(window_size))).multiply(impulse_response).reduce("sum", [0]).get([0]);
    });
    return smooth;
  },
  
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
  
  addMultiSmoothImageDoySeries: function(roi, dataset, bandNames, pixelScale) {
    var raw_data_list = dataset.select(bandNames).getRegion(roi, pixelScale);
    var bandCount = bandNames.length;
    raw_data_list.evaluate(function(result) {
      var bands_value_list = [];
      for(var j=0; j<bandCount; j++) {
        bands_value_list.push([]);
      }
      
      var doy_list = [];
      for (var i =0; i<result.length; i++) {
        var data = result[i];
        if (data[0] !== "id") {
          var time_str = ee.Date(data[3]).format("DDD");
          doy_list.push(time_str);
          for(var k=0; k<bandCount; k++) {
            bands_value_list[k].push(data[4+k]);
          }
        }
      }
      
      var smooth_data_list = [];
      for(j=0; j<bandCount; j++) {
        smooth_data_list.push([]);
      }
      for (i = 0; i < bands_value_list.length; i++) {
        var d_list = bands_value_list[i];
        for (j=0; j<d_list.length; j++) {
          var d = d_list[j];
          if (d === null) {
            d = SG_Models.fill_nodata(j, d_list, 3);
          }
          smooth_data_list[i].push(d);
        }
      }

      var sg_filter_list = [];
      var window_size = 7;
      var order = 2;
      for (i=0; i<smooth_data_list.length; i++) {
        var src_list = ee.List(smooth_data_list[i]);
        var smooth = SG_Models.sgFilter(src_list, window_size, order);
        sg_filter_list.push(smooth);
      }
      var yValues = ee.Array.cat(sg_filter_list, 1);
      var chart = ui.Chart.array.values(yValues, 0, ee.List(doy_list))
                    .setSeriesNames(bandNames)
                    .setOptions({
                      title: 'Window: ' + window_size + ' Order: ' + order, 
                      hAxis: {title: 'time'}, vAxis: {title: "Value"},
                      legend: null,
                      lineWidth:1,
                      pointSize:2
                    });
      print(chart);
    });
  }
};



exports.getS2Collection = function(startDate, endDate, studyArea, bandNames) {
  var dataset = S2.filterS2Data(startDate, endDate, studyArea)
                  .select(bandNames);
  return dataset;
};

exports.getS2RawCollection = function(startDate, endDate, studyArea) {
  var dataset = S2.filterS2RawData(startDate, endDate, studyArea);
  return dataset;
};

exports.s2RCLE = function(img) {
  return S2.RCLE(img);
};

///////////////////////////////////////////////////////////
exports.addMultiSmoothImageDoySeries = function(roi, dataset, bandNames, pixelScale) {
  SG_Models.addMultiSmoothImageDoySeries(roi, dataset, bandNames, pixelScale);
};


