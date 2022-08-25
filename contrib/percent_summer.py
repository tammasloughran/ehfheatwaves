#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 16:32:54 2018

@author: Tammas Loughran
"""
import numpy as np
import netCDF4 as nc
import pandas as pd
import datetime as dt
ncfile = nc.Dataset('HadGHCND_TXTN_1950-2014.nc','r')
data = ncfile.variables['tmax'][:]
if type(data)!=np.ma.core.MaskedArray: data = np.ma.array(data, mask=np.isnan(data))
data = np.logical_not(data.mask)
time = ncfile.variables['time'][:]
dayone = dt.datetime(int(str(time[0])[:4]),int(str(time[0])[4:6]),int(str(time[0])[6:8]))
daylast = dt.datetime(int(str(time[-1])[:4]),int(str(time[-1])[4:6]),int(str(time[-1])[6:8]))
dates = pd.date_range(dayone,daylast)
percent = np.ones((daylast.year-dayone.year+1, data.shape[1], data.shape[2]))*np.nan
for i,year in enumerate(range(dayone.year,daylast.year+1)):
    year_index = (dates.year==year)&((dates.month==11)|(dates.month==12)|(dates.month==1)|(dates.month==2)|(dates.month==13))
    percent[i,...] = 100*data[year_index,...].sum(axis=0)/year_index.sum()
outfile = nc.Dataset('percent.nc','w')
outfile.createDimension('time', size=percent.shape[0])
outfile.createDimension('latitude',size=ncfile.dimensions['latitude'].size)
outfile.createDimension('longitude',size=ncfile.dimensions['longitude'].size)
timevar = outfile.createVariable('time', float, dimensions=('time'))
latsvar = outfile.createVariable('latitude', float, dimensions=('latitude'))
lonsvar = outfile.createVariable('longitude', float, dimensions=('longitude'))
percentvar = outfile.createVariable('percent', float, dimensions=('time','latitude','longitude'))
setattr(timevar, 'units', 'years since 1949-01-01')
setattr(latsvar, 'standard_name', 'latitude')
setattr(latsvar, 'units', 'degrees_north')
setattr(latsvar, 'axis', 'Y')
setattr(lonsvar, 'standard_name', 'longitude')
setattr(lonsvar, 'units', 'degrees_east')
setattr(lonsvar, 'axis', 'X')
timevar[:] = np.array(range(dayone.year, daylast.year+1)) - dayone.year
latsvar[:] = ncfile.variables['latitude'][:]
lonsvar[:] = ncfile.variables['longitude'][:]
percentvar[:] = percent
outfile.close()