# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 07:20:28 2020

@author: Administrator
"""

import ee
from ipygee import *
import pandas as pd
import numpy as np
import os
from datetime import datetime

site_pd = pd.read_csv(r"***\loc_info.csv")
outPath = r"***" 
ee.Initialize()
pixelScale = 1000 #buffer
CFSV2 = ee.ImageCollection("NOAA/CFSV2/FOR6H")

startDate = "2018-1-1"
endDate = "2019-1-1"

dataset = (CFSV2.filterDate(startDate, endDate)
               .sort("system:time_start")) 

band = ['Precipitation_rate_surface_6_Hour_Average']

for s in range(len(site_pd)):
    site_name = site_pd.iloc[s,0] #siteName
    print(site_name)
    lon = site_pd.iloc[s,1]
    lat = site_pd.iloc[s,2]
    
    point = ee.Geometry.Point(lon, lat)
    
    ts = chart.Image.series(**{
        'imageCollection': dataset.select(band),  
        'region': point,
        'reducer': ee.Reducer.mean(),
        'scale': pixelScale,
        'bands': band,
        'label_bands':['mean_Pre'],
    })
    
    data = ts.dataframe
    time_step = 4
    #6-Hour Data 
    arr = np.full((int(len(data)/time_step), 2), -9999.0)    

    for i in range(int(len(data)/time_step)):
        total = 0
        for j in range(0, time_step):
            n = i*time_step + j
            total = total + data.iloc[n,0]            
        arr[i,0] = i + 1
        arr[i,1] = (total/time_step)*6*60*60   
    
    outfile = os.path.join(outPath, "{0}.csv".format(site_name))
    np.savetxt(outfile, arr, delimiter = ',')
    