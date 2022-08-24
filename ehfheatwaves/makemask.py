#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Tammas Loughran
"""
import netCDF4 as nc
import numpy as np
# Load file
ncfile = nc.Dataset('HadGHCND_TXTN_1950-2014.nc','r')
data = ncfile.variables['tmax'][:]
# Make mask where data exists at any point in time
mask = np.array(data.any(axis=0).astype(int))
# Save to mask.nc
maskfile = nc.Dataset('mask.nc','w')
maskfile.createDimension('latitude', size=mask.shape[0])
maskfile.createDimension('longitude', size=mask.shape[1])
latsvar = maskfile.createVariable('latitude', float, dimensions=('latitude'))
lonsvar = maskfile.createVariable('longitude', float, dimensions=('longitude'))
maskvar = maskfile.createVariable('mask', int, dimensions=('latitude','longitude'))
setattr(latsvar, 'standard_name', 'latitude')
setattr(latsvar, 'units', 'degrees_north')
setattr(latsvar, 'axis', 'Y')
setattr(lonsvar, 'standard_name', 'longitude')
setattr(lonsvar, 'units', 'degrees_east')
setattr(lonsvar, 'axis', 'X')
latsvar[:] = ncfile.variables['latitude'][:]
lonsvar[:] = ncfile.variables['longitude'][:]
maskvar[:] = mask
ncfile.close()
maskfile.close()