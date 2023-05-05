#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
ncio.py contains the functions for reading from and writing to NetCDF4 files.
Created on Sat Apr 14 12:56:42 2018

@author: Tammas Loughran
"""
import datetime as dt
import sys

import netCDF4 as nc
import numpy as np
import pandas as pd

import ehfheatwaves.constants as const
from ehfheatwaves.__init__ import __version__

LAT_NAMES = ('lat','lats','latitude','latitudes')
LON_NAMES = ('lon','lons','longitude','longitudes')


class DatesOrderError(Exception):
    """Exception to be raised when the calendar start day occurs before the end day."""

    def __init__(self, start, end):
        self.start = start
        self.end = end
        print("Calendar end date appears before the start date")
        print("start: ", self.start)
        print("end: ", self.end)


class Calendar360():
    """A basic calendar class with 360 days per year."""

    def __init__(self, sdate, edate):
        """## Arguments:

        - sdate : Start date of calendar.
        - edate : End date of calendar.
        """
        if sdate>edate: raise DatesOrderError(sdate, edate)
        self.year = np.repeat(range(sdate.year, edate.year+1), 360, 0)
        nyears = len(range(sdate.year, edate.year+1))
        self.month = np.tile(np.repeat(range(1, 12+1), 30, 0), nyears)
        self.day = np.tile(np.tile(range(1, 30+1), 12), nyears)
        if (sdate.day!=1)|(edate.month!=1):
            sdoyi = (sdate.month - 1)*30 + sdate.day - 1
            self.year = self.year[sdoyi:]
            self.month = self.month[sdoyi:]
            self.day = self.day[sdoyi:]
        if (edate.day!=30)|(edate.month!=12):
            edoyi = (12 - edate.month)*30 + (30 - edate.day)
            self.year = self.year[:-edoyi]
            self.month = self.month[:-edoyi]
            self.day = self.day[:-edoyi]

    def __getitem__(self, i):
        """Return a datetime object for the provided index."""
        return dt.datetime(self.year[i], self.month[i], self.day[i])


class TimeData(object):
    """Data class for time data from an netcdf file."""


    def __init__(self, filename:str, timevname:str):
        """## Arguments:

        - filename : Input file name.
        - timevarname : Name of the time variable in the input file.
        """
        if any([(wildcard in filename) for wildcard in ['*','?','[']]):
            tempnc = nc.MFDataset(filename, 'r')
            nctime = tempnc.variables[timevname]
            nctime = nc.MFTime(nctime)
        else:
            tempnc = nc.Dataset(filename, 'r')
            nctime = tempnc.variables[timevname]

        try:
            self.calendar = nctime.calendar
        except:
            self.calendar = 'proleptic_gregorian'
        if not self.calendar:
            print('Unrecognized calendar. Using gregorian.')
            self.calendar = 'gregorian'

        if self.calendar=='360_day':
            self.daysinyear = 360
            # 360 day season start and end indices
            self.SHS = (300,450)
            self.SHW = (120,270)
            self.dayone = nc.num2date(nctime[0], nctime.units, calendar=self.calendar)
            self.daylast = nc.num2date(nctime[-1], nctime.units, calendar=self.calendar)
            self.dates = Calendar360(self.dayone, self.daylast)
        else:
            self.daysinyear = 365
            # 365 day season start and end indices
            self.SHS = (304,455)
            self.SHW = (120,273)
            if tempnc.variables[timevname].units=='day as %Y%m%d.%f':
                st = str(int(nctime[0]))
                nd = str(int(nctime[-1]))
                self.dayone = dt.datetime(int(st[:4]), int(st[4:6]), int(st[6:]))
                self.daylast = dt.datetime(int(nd[:4]), int(nd[4:6]), int(nd[6:]))
            else:
                self.dayone = nc.num2date(nctime[0], nctime.units, calendar=self.calendar)
                self.daylast = nc.num2date(nctime[-1], nctime.units, calendar=self.calendar)
            self.dates = pd.period_range(str(self.dayone), str(self.daylast))
            # Remove leap days. Maybe this should be a separate function?
            self.noleapdates = self.dates[(self.dates.month!=2)|(self.dates.day!=29)]


def get_mask(options)->np.ndarray:
    """get_mask loads the land-sea mask data from a NetCDF4 file."""
    masknc = nc.Dataset(options.maskfile, 'r')
    mask = masknc.variables[options.maskvname][:]
    if mask.max()>1: mask = mask>50
    else: mask = mask>0.5
    mask = mask.astype(bool)
    mask = np.squeeze(mask)
    masknc.close()
    if options.invertmask: mask = np.logical_not(mask)
    if options.flipmask: mask = np.flipud(mask)
    return mask


def load_bp_data(options, timedata:TimeData, variable:str='tmax', mask:np.ndarray=None)->np.ndarray:
    """Load the tmax or tmin data for the base period.

    ## Arguments:

    - options
    - timedata : `TimeData` object.
    - variable : The variable to load.
    - mask : The land-sea mask.

    ## Returns:

    - temp : Temperature data array.
    """
    # Determine if we need tmax or tmin and whether the data are in multiple files.
    if variable=='tmax':
        if options.bpfx:files = options.bpfx
        else: files = options.tmaxfile
        varname = options.tmaxvname
    elif variable=='tmin':
        if options.bpfn: files = options.bpfn
        else: files = options.tminfile
        varname = options.tminvname

    # Load ncfile depending on multi or single-file
    if any([(wildcard in files) for wildcard in ['*','?','[']]):
        tempnc = nc.MFDataset(files, 'r')
        bptime = tempnc.variables[options.timevname]
        bptime = MFTime(bptime)
    else:
        tempnc = nc.Dataset(files,'r')
        bptime = tempnc.variables[options.timevname]

    # Define the base period
    if tempnc.variables[options.timevname].units=='day as %Y%m%d.%f':
        st = str(int(bptime[0]))
        nd = str(int(bptime[-1]))
        bpdayone = dt.datetime(int(st[:4]), int(st[4:6]), int(st[6:]))
        bpdaylast = dt.datetime(int(nd[:4]), int(nd[4:6]), int(nd[6:]))
    else:
        bpdayone = nc.num2date(bptime[0], bptime.units, calendar=timedata.calendar)
        bpdaylast = nc.num2date(bptime[-1], bptime.units, calendar=timedata.calendar)

    if timedata.calendar=='360_day': bpdates = Calendar360(bpdayone, bpdaylast)
    else:
        bpdates = pd.period_range(str(bpdayone), str(bpdaylast))
        if timedata.calendar=='365_day': bpdates = bpdates[(bpdates.month!=2)|(bpdates.day!=29)]
        dates_base = bpdates[(options.bpstart<=bpdates.year)&(bpdates.year<=options.bpend)]

    try:
        temp = tempnc.variables[varname][(options.bpstart<=bpdates.year)&(bpdates.year<=options.bpend)]
    except IndexError as FirstException:
        raise IndexError(
                "Could not index netCDF file time dimension. This could be a calendar problem. "
                "Please make sure your input file is daily frequency data and/or that it contains "
                f"the base period. You wanted {options.bpstart}-{options.bpend} and your data have "
                f"{bpdates.year[0]}-{bpdates.year[-1]}"
                ) from FirstException
    if len(temp.shape)==4: temp = temp.squeeze()

    if options.maskfile:
        temp = temp[:,mask]

    if tempnc.variables[varname].units=='K': temp -= 273.15

    # Remove leap days in gregorian calendars
    if (timedata.calendar=='gregorian')|(timedata.calendar=='proleptic_gregorian')|(timedata.calendar=='standard'):
        temp = temp[(dates_base.month!=2)|(dates_base.day!=29),...]

    tempnc.close()
    return temp


def remove_leap_days(data:np.ndarray, dates:pd.DatetimeIndex)->np.ndarray:
    """Remove February 29th from a dataset."""
    return data[(dates.month!=2)|(dates.day!=29),...]


def get_all_data(files:str, vname:str, options)->tuple:
    """Load all temperature data from a netcdf file.

    ## Arguments:

    - files : filename, may contain wildcards.
    - vname : variable name.
    - options

    ## Returns:

    - temp : data in (time, x, y) coordinates.
    - lats : latitudes
    """
    if any([(wildcard in files) for wildcard in ['*','?','[']]):
        tempnc = nc.MFDataset(files, 'r')
    else:
        tempnc = nc.Dataset(files, 'r')
    temp = tempnc.variables[vname][:]
    if len(temp.shape)==4: temp = temp.squeeze()
    # Test for increasing latitude and flip if decreasing
    latkey = [vrbl in LAT_NAMES for vrbl in tempnc.variables.keys()].index(True)
    latvname = list(tempnc.variables.keys())[latkey]
    lats = tempnc.variables[latvname][:]
    if (lats[-1]-lats[0])<0: lats = np.flipud(lats)
    if tempnc.variables[vname].units=='K': temp -= 273.15
    tempnc.close()
    return temp, lats


def save_yearly(
            HWA,
            HWM,
            HWN,
            HWF,
            HWD,
            HWT,
            tpct,
            definition:str,
            timedata,
            options,
            mask:np.ndarray,
            )->None:
    """Save yearly data to netcdf file.

    Input aspect arrays are 2D, timeXspace and are either reshaped or indexed with a land-sea mask.

    ## Arguments:

    - HWA :
    - HWM :
    - HWN :
    - HWF :
    - HWD :
    - HWT :
    - txpct : Percentile thresholds.
    - definition : Name of heatwave definition.
    - timedata : TimeData onbject.
    - options :
    - mask : Land-sea mask.
    """
    if any([(wildcard in options.tmaxfile) for wildcard in ['*','?','[']]):
        tempnc = nc.MFDataset(options.tmaxfile, 'r')
    else:
        tempnc = nc.Dataset(options.tmaxfile, 'r')
    nyears = HWA.shape[0]
    try: experiment = tempnc.__getattribute__('experiment')
    except AttributeError: experiment = ''
    try: model = tempnc.__getattribute__('model_id')
    except AttributeError: model = ''
    try: parent = tempnc.__getattribute__('parent_experiment_rip')
    except AttributeError: parent = ''
    try: realization = tempnc.__getattribute__('realization')
    except AttributeError: realization = ''
    try: initialization = tempnc.__getattribute__('initialization_method')
    except AttributeError: initialization = ''
    try:
        physics = tempnc.__getattribute__('physics_version')
        rip = 'r'+str(realization)+'i'+str(initialization)+'p'+str(physics)
    except AttributeError:
        physics = ''
        rip = ''
    latkey = [vrbl in LAT_NAMES for vrbl in tempnc.variables.keys()].index(True)
    latvname = list(tempnc.variables.keys())[latkey]
    lonkey = [vrbl in LON_NAMES for vrbl in tempnc.variables.keys()].index(True)
    lonvname = list(tempnc.variables.keys())[lonkey]
    space = (tempnc.dimensions[latvname].__len__(), tempnc.dimensions[lonvname].__len__())
    yearlyout = nc.Dataset('%s_heatwaves_%s_%s_%s_yearly_%s.nc'%(definition, model, experiment, rip, options.season), 'w')
    yearlyout.createDimension('time', size=None)
    yearlyout.createDimension('lon', tempnc.dimensions[lonvname].__len__())
    yearlyout.createDimension('lat', tempnc.dimensions[latvname].__len__())
    yearlyout.createDimension('day', timedata.daysinyear)
    setattr(yearlyout, "author", "Tammas Loughran")
    setattr(yearlyout, "contact", "t.loughran@unsw.edu.au")
    setattr(yearlyout, "source", "https://github.com/tammasloughran/ehfheatwaves")
    setattr(yearlyout, "date", dt.datetime.today().strftime('%Y-%m-%d'))
    setattr(yearlyout, "script", sys.argv[0])
    if model:
        setattr(yearlyout, "model_id", model)
        setattr(yearlyout, "experiment", experiment)
        setattr(yearlyout, "parent_experiment_rip", parent)
        setattr(yearlyout, "realization", realization)
        setattr(yearlyout, "initialization_method", initialization)
        setattr(yearlyout, "physics_version", physics)
    setattr(yearlyout, "period", "%s-%s"%(str(timedata.dayone.year),str(timedata.daylast.year)))
    setattr(yearlyout, "base_period", "%s-%s"%(str(options.bpstart),str(options.bpend)))
    setattr(yearlyout, "percentile", "%sth"%(str(options.pcntl)))
    setattr(yearlyout, "definition", definition)
    setattr(yearlyout, "frequency", "yearly")
    setattr(yearlyout, "season", options.season)
    setattr(yearlyout, "season_note", ("The year of a season is the year it starts in. SH summer: Nov-Mar. NH summer: May-Sep."))
    setattr(yearlyout, "version", __version__)
    setattr(yearlyout, "tmax_file", options.tmaxfile)
    setattr(yearlyout, "tmin_file", options.tminfile)
    if options.maskfile:
        setattr(yearlyout, "mask_file", options.maskfile)
    setattr(yearlyout, "quantile_method", options.qtilemethod)
    setattr(yearlyout, 'Conventions', 'CF-1.7')
    otime = yearlyout.createVariable('time', 'f8', 'time')
    setattr(otime, 'units', 'years since 0001-06-01 00:00:00')
    setattr(otime, 'standard_name', 'time')
    setattr(otime, 'calendar', 'proleptic_gregorian')
    olat = yearlyout.createVariable('lat', 'f8', 'lat')
    setattr(olat, 'standard_name', 'latitude')
    setattr(olat, 'long_name', 'Latitude')
    setattr(olat, 'units', 'degrees_north')
    setattr(olat, 'axis', 'Y')
    olon = yearlyout.createVariable('lon', 'f8', 'lon')
    setattr(olon, 'standard_name', 'longitude')
    setattr(olon, 'long_name', 'Longitude')
    setattr(olon, 'units', 'degrees_east')
    setattr(olon, 'axis', 'X')
    otpct = yearlyout.createVariable('t%spct'%(options.pcntl), 'f8', ('day','lat','lon'), fill_value=const.FILL_VAL)
    setattr(otpct, 'long_name', '90th percentile')
    setattr(otpct, 'units', 'degC')
    setattr(otpct, 'description', '90th percentile of %s-%s'%(str(options.bpstart),str(options.bpend)))
    setattr(otpct, 'missing_value', const.MISSING_VAL)
    setattr(otpct, 'valid_range', (-20,100))
    HWAout = yearlyout.createVariable('HWA_%s'%(definition), 'f8', ('time','lat','lon'), fill_value=const.FILL_VAL)
    setattr(HWAout, 'long_name', 'Heatwave Amplitude')
    if definition=='tx90pct':
        setattr(HWAout, 'units', 'degC')
    elif definition=='tn90pct':
        setattr(HWAout, 'units', 'degC')
    elif definition=='EHF':
        setattr(HWAout, 'units', 'degC2')
    setattr(HWAout, 'description', 'Peak of the hottest heatwave per year')
    setattr(HWAout, 'missing_value', const.MISSING_VAL)
    setattr(HWAout, 'valid_range', (-30, 400))
    HWMout = yearlyout.createVariable('HWM_%s'%(definition), 'f8', ('time','lat','lon'), fill_value=const.FILL_VAL)
    setattr(HWMout, 'long_name', 'Heatwave Magnitude')
    if definition=='tx90pct':
        setattr(HWMout, 'units', 'degC')
    elif definition=='tn90pct':
        setattr(HWMout, 'units', 'degC')
    elif definition=='EHF':
        setattr(HWMout, 'units', 'degC2')
    setattr(HWMout, 'description', 'Average magnitude of the yearly heatwave')
    setattr(HWMout, 'missing_value', const.MISSING_VAL)
    setattr(HWMout, 'valid_range' ,(-30, 400))
    HWNout = yearlyout.createVariable('HWN_%s'%(definition), 'f8', ('time', 'lat', 'lon'), fill_value=const.FILL_VAL)
    setattr(HWNout, 'long_name', 'Heatwave Number')
    setattr(HWNout, 'units','count')
    setattr(HWNout, 'description', 'Number of heatwaves per year')
    setattr(HWNout, 'missing_value', const.MISSING_VAL)
    setattr(HWNout, 'valid_range', (0, 40))
    HWFout = yearlyout.createVariable('HWF_%s'%(definition), 'f8', ('time','lat','lon'), fill_value=const.FILL_VAL)
    setattr(HWFout, 'long_name','Heatwave Frequency')
    setattr(HWFout, 'units', 'days')
    setattr(HWFout, 'description', 'Proportion of heatwave days per season')
    setattr(HWFout, 'missing_value', const.MISSING_VAL)
    setattr(HWFout, 'valid_range', (0, 165))
    HWDout = yearlyout.createVariable('HWD_%s'%(definition), 'f8', ('time','lat','lon'), fill_value=const.FILL_VAL)
    setattr(HWDout, 'long_name', 'Heatwave Duration')
    setattr(HWDout, 'units', 'days')
    setattr(HWDout, 'description', 'Duration of the longest heatwave per year')
    setattr(HWDout, 'missing_value', const.MISSING_VAL)
    setattr(HWDout, 'valid_range', (0,165))
    HWTout = yearlyout.createVariable('HWT_%s'%(definition), 'f8', ('time','lat','lon'), fill_value=const.FILL_VAL)
    setattr(HWTout, 'long_name', 'Heatwave Timing')
    if options.season=='summer':
        setattr(HWTout, 'units', 'days since 0001-11-01 00:00:00')
    elif options.season=='winter':
        setattr(HWTout, 'units', 'days since 0001-05-01 00:00:00')
    setattr(HWTout, 'description', 'First heat wave day of the season')
    setattr(HWTout, 'missing_value', const.MISSING_VAL)
    setattr(HWTout, 'valid_range', (1,366))
    otime[:] = range(timedata.dayone.year, timedata.daylast.year)
    olat[:] = tempnc.variables[latvname][:]
    olon[:] = tempnc.variables[lonvname][:]
    if options.maskfile:
        dummy_array = np.ones((timedata.daysinyear,)+(len(olat),)+(len(olon),))*const.FILL_VAL
        dummy_array[:,mask] = tpct
        otpct[:] = dummy_array.copy()
        dummy_array = np.ones((nyears,)+(len(olat),)+(len(olon),))*const.FILL_VAL
        dummy_array[:,mask] = HWA
        HWAout[:] = dummy_array.copy()
        dummy_array = np.ones((nyears,)+(len(olat),)+(len(olon),))*const.FILL_VAL
        dummy_array[:,mask] = HWM
        HWMout[:] = dummy_array.copy()
        dummy_array = np.ones((nyears,)+(len(olat),)+(len(olon),))*const.FILL_VAL
        dummy_array[:,mask] = HWN
        HWNout[:] = dummy_array.copy()
        dummy_array = np.ones((nyears,)+(len(olat),)+(len(olon),))*const.FILL_VAL
        dummy_array[:,mask] = HWF
        HWFout[:] = dummy_array.copy()
        dummy_array = np.ones((nyears,)+(len(olat),)+(len(olon),))*const.FILL_VAL
        dummy_array[:,mask] = HWD
        HWDout[:] = dummy_array.copy()
        dummy_array = np.ones((nyears,)+(len(olat),)+(len(olon),))*const.FILL_VAL
        dummy_array[:,mask] = HWT
        HWTout[:] = dummy_array.copy()
    else:
        otpct[:] = tpct
        HWAout[:] = HWA.reshape((nyears,)+space)
        HWMout[:] = HWM.reshape((nyears,)+space)
        HWNout[:] = HWN.reshape((nyears,)+space)
        HWFout[:] = HWF.reshape((nyears,)+space)
        HWDout[:] = HWD.reshape((nyears,)+space)
        HWTout[:] = HWT.reshape((nyears,)+space)
    tempnc.close()
    yearlyout.close()


def save_daily(
            exceed:np.ndarray,
            event:np.ndarray,
            ends:np.ndarray,
            options,
            timedata:TimeData,
            original_shape:tuple,
            mask:np.ndarray,
            defn:str='EHF',
            )->None:
    """save_daily saves the daily data to netcdf file. Input arrays are 2D, (time, space) and are
    either reshaped or indexed with a land-sea mask.

    ## Arguments:

    - exceed : Threshold exceedence values.
    - event : Boolean indicator for whether a heatwave is occuring.
    - ends : Duration of heatwave.
    - options :
    - timedata : TimeData object.
    - original_shape : The original shape of the input data.
    - mask : Land-sea mask.
    - defn : The name of the heatwave definition being saved.
    """
    if options.tmaxfile:
        filename = options.tmaxfile
    elif options.tminfile:
        filename = options.tminfile
    if any([(wildcard in filename) for wildcard in ['*','?','[']]):
        tempnc = nc.MFDataset(filename, 'r')
    else:
        tempnc = nc.Dataset(filename, 'r')
    try: experiment = tempnc.__getattribute__('experiment')
    except AttributeError: experiment = ''
    try: model = tempnc.__getattribute__('model_id')
    except AttributeError: model = ''
    try: parent = tempnc.__getattribute__('parent_experiment_rip')
    except AttributeError: parent = ''
    try: realization = tempnc.__getattribute__('realization')
    except AttributeError: realization = ''
    try: initialization = tempnc.__getattribute__('initialization_method')
    except AttributeError: initialization = ''
    try:
        physics = tempnc.__getattribute__('physics_version')
        rip = 'r'+str(realization)+'i'+str(initialization)+'p'+str(physics)
    except AttributeError:
        physics = ''
        rip = ''
    latkey = [vrbl in LAT_NAMES for vrbl in tempnc.variables.keys()].index(True)
    latvname = list(tempnc.variables.keys())[latkey]
    lonkey = [vrbl in LON_NAMES for vrbl in tempnc.variables.keys()].index(True)
    lonvname = list(tempnc.variables.keys())[lonkey]
    dailyout = nc.Dataset('%s_heatwaves_%s_%s_%s_daily.nc'%(defn, model, experiment, rip), mode='w')
    dailyout.createDimension('time', size=None)
    dailyout.createDimension('lon', tempnc.dimensions[lonvname].__len__())
    dailyout.createDimension('lat', tempnc.dimensions[latvname].__len__())
    setattr(dailyout, "author", "Tammas Loughran")
    setattr(dailyout, "contact", "t.loughran@unsw.edu.au")
    setattr(dailyout, "source", "https://github.com/tammasloughran/ehfheatwaves")
    setattr(dailyout, "date", dt.datetime.today().strftime('%Y-%m-%d'))
    setattr(dailyout, "script", sys.argv[0])
    setattr(dailyout, "period", "%s-%s"%(str(timedata.dayone.year), str(timedata.daylast.year)))
    setattr(dailyout, "base_period", "%s-%s"%(str(options.bpstart), str(options.bpend)))
    setattr(dailyout, "percentile", "%sth"%(str(options.pcntl)))
    if model:
        setattr(dailyout, "model_id", model)
        setattr(dailyout, "experiment", experiment)
        setattr(dailyout, "parent_experiment_rip", parent)
        setattr(dailyout, "realization", realization)
        setattr(dailyout, "initialization_method", initialization)
        setattr(dailyout, "physics_version", physics)
    setattr(dailyout, "version", __version__)
    if options.tmaxfile:
        setattr(dailyout, "tmax_file", options.tmaxfile)
    if options.tminfile:
        setattr(dailyout, "tmin_file", options.tminfile)
    if options.maskfile:
        setattr(dailyout, "mask_file", str(options.maskfile))
    setattr(dailyout, "quantile_method", options.qtilemethod)
    setattr(dailyout, 'Conventions', 'CF-1.7')
    otime = dailyout.createVariable('time', 'f8', 'time')
    setattr(otime, 'units', 'days since %s-01-01'%(timedata.dayone.year))
    setattr(otime, 'calendar', '365_day')
    setattr(otime, 'standard_name', 'time')
    olat = dailyout.createVariable('lat', 'f8', 'lat')
    setattr(olat, 'standard_name', 'latitude')
    setattr(olat, 'long_name', 'Latitude')
    setattr(olat, 'units', 'degrees_north')
    olon = dailyout.createVariable('lon', 'f8', 'lon')
    setattr(olon, 'standard_name', 'longitude')
    setattr(olon, 'long_name', 'Longitude')
    setattr(olon, 'units', 'degrees_east')
    oehf = dailyout.createVariable(defn, 'f8', ('time','lat','lon'), fill_value=const.FILL_VAL)
    if defn=='EHF':
        setattr(oehf, 'long_name', 'Excess Heat Factor')
        setattr(oehf, 'units', 'degC2')
    elif defn=='tx90pct':
        setattr(oehf, 'long_name', 'Temperature Exceeding tx90pct')
        setattr(oehf, 'units', 'C')
    elif defn=='tn90pct':
        setattr(oehf, 'long_name', 'Temperature Exceeding tn90pct')
        setattr(oehf, 'units', 'C')
    setattr(oehf, 'missing_value', const.MISSING_VAL)
    oevent = dailyout.createVariable('event', 'f8', ('time','lat','lon'), fill_value=const.FILL_VAL)
    setattr(oevent, 'long_name', 'Event indicator')
    setattr(oevent, 'description', 'Indicates whether a heatwave is happening on that day')
    setattr(oevent, 'missing_value', const.MISSING_VAL)
    oends = dailyout.createVariable('ends', 'f8', ('time','lat','lon'), fill_value=const.FILL_VAL)
    setattr(oends, 'long_name', 'Duration at start of heatwave')
    setattr(oends, 'units', 'days')
    setattr(oends, 'missing_value', const.MISSING_VAL)
    otime[:] = range(0,original_shape[0],1)
    olat[:] = tempnc.variables[latvname][:]
    olon[:] = tempnc.variables[lonvname][:]
    exceed[exceed.mask==True] = const.MISSING_VAL
    if options.maskfile:
        dummy_array = np.ones(original_shape)*const.FILL_VAL
        dummy_array[:,mask] = exceed
        oehf[:] = dummy_array.copy()
        dummy_array = np.ones(original_shape)*const.FILL_VAL
        dummy_array[:,mask] = event
        oevent[:] = dummy_array.copy()
        dummy_array = np.ones(original_shape)*const.FILL_VAL
        dummy_array[:,mask] = ends
        oends[:] = dummy_array.copy()
    else:
        oehf[:] = exceed
        oevent[:] = event
        oends[:] = ends
    tempnc.close()
    dailyout.close()


def save_ehi(
            EHIsig:np.ndarray,
            EHIaccl:np.ndarray,
            options,
            timedata:TimeData,
            original_shape:tuple,
            mask:np.ndarray,
            )->None:
    """Save the daily data to netcdf file. Input arrays are 2D, (time, space) and are either
    reshaped or indexed with a land-sea mask.

    ## Arguments:

    - EHIsig : EHI significance index.
    - EHIaccl : EHI acclimatisation index.
    - options :
    - timedata : TimeData object.
    - original_shape : The original shape of the input data.
    - mask : Land-sea mask.
    """
    if options.tmaxfile:
        filename = options.tmaxfile
    elif options.tminfile:
        filename = options.tminfile
    if any([(wildcard in filename) for wildcard in ['*','?','[']]):
        tempnc = nc.MFDataset(filename, 'r')
    else:
        tempnc = nc.Dataset(filename, 'r')
    try: experiment = tempnc.__getattribute__('experiment')
    except AttributeError: experiment = ''
    try: model = tempnc.__getattribute__('model_id')
    except AttributeError: model = ''
    try: parent = tempnc.__getattribute__('parent_experiment_rip')
    except AttributeError: parent = ''
    try: realization = tempnc.__getattribute__('realization')
    except AttributeError: realization = ''
    try: initialization = tempnc.__getattribute__('initialization_method')
    except AttributeError: initialization = ''
    try:
        physics = tempnc.__getattribute__('physics_version')
        rip = 'r'+str(realization)+'i'+str(initialization)+'p'+str(physics)
    except AttributeError:
        physics = ''
        rip = ''
    latkey = [vrbl in LAT_NAMES for vrbl in tempnc.variables.keys()].index(True)
    latvname = list(tempnc.variables.keys())[latkey]
    lonkey = [vrbl in LON_NAMES for vrbl in tempnc.variables.keys()].index(True)
    lonvname = list(tempnc.variables.keys())[lonkey]
    defn = 'EHI'
    dailyout = nc.Dataset('%s_heatwaves_%s_%s_%s_daily.nc'%(defn, model, experiment, rip), mode='w')
    dailyout.createDimension('time', size=None)
    dailyout.createDimension('lon', tempnc.dimensions[lonvname].__len__())
    dailyout.createDimension('lat', tempnc.dimensions[latvname].__len__())
    setattr(dailyout, "author", "Tammas Loughran")
    setattr(dailyout, "contact", "t.loughran@unsw.edu.au")
    setattr(dailyout, "source", "https://github.com/tammasloughran/ehfheatwaves")
    setattr(dailyout, "date", dt.datetime.today().strftime('%Y-%m-%d'))
    setattr(dailyout, "script", sys.argv[0])
    setattr(dailyout, "period", "%s-%s"%(str(timedata.dayone.year), str(timedata.daylast.year)))
    setattr(dailyout, "base_period", "%s-%s"%(str(options.bpstart), str(options.bpend)))
    setattr(dailyout, "percentile", "%sth"%(str(options.pcntl)))
    if model:
        setattr(dailyout, "model_id", model)
        setattr(dailyout, "experiment", experiment)
        setattr(dailyout, "parent_experiment_rip", parent)
        setattr(dailyout, "realization", realization)
        setattr(dailyout, "initialization_method", initialization)
        setattr(dailyout, "physics_version", physics)
    setattr(dailyout, "version", __version__)
    if options.tmaxfile:
        setattr(dailyout, "tmax_file", options.tmaxfile)
    if options.tminfile:
        setattr(dailyout, "tmin_file", options.tminfile)
    if options.maskfile:
        setattr(dailyout, "mask_file", str(options.maskfile))
    setattr(dailyout, "quantile_method", options.qtilemethod)
    setattr(dailyout, 'Conventions', 'CF-1.7')
    otime = dailyout.createVariable('time', 'f8', 'time')
    setattr(otime, 'standard_name', 'time')
    setattr(otime, 'units', 'days since %s-01-01'%(timedata.dayone.year))
    setattr(otime, 'calendar', '365_day')
    olat = dailyout.createVariable('lat', 'f8', 'lat')
    setattr(olat, 'standard_name', 'latitude')
    setattr(olat, 'long_name', 'Latitude')
    setattr(olat, 'units', 'degrees_north')
    olon = dailyout.createVariable('lon', 'f8', 'lon')
    setattr(olon, 'standard_name', 'longitude')
    setattr(olon, 'long_name', 'Longitude')
    setattr(olon, 'units', 'degrees_east')
    oehisig = dailyout.createVariable(
            'EHIsig',
            'f8',
            ('time','lat','lon'),
            fill_value=const.FILL_VAL,
            )
    setattr(oehisig, 'long_name', 'Excess Heat Index Significance')
    setattr(oehisig, 'units', 'C')
    setattr(oehisig, 'missing_value', const.MISSING_VAL)
    setattr(oehisig, 'valid_range', (0,100))
    oehiaccl = dailyout.createVariable(
            'EHIaccl',
            'f8',
            ('time','lat','lon'),
            fill_value=const.FILL_VAL,
            )
    setattr(oehiaccl, 'long_name', 'Excess Heat Index Acclimatisation')
    setattr(oehiaccl, 'units', 'C')
    setattr(oehiaccl, 'missing_value', const.MISSING_VAL)
    setattr(oehiaccl, 'valid_range', (0,100))
    otime[:] = range(0, original_shape[0], 1)
    olat[:] = tempnc.variables[latvname][:]
    olon[:] = tempnc.variables[lonvname][:]
    if options.maskfile:
        dummy_array = np.ones(original_shape)*const.FILL_VAL
        dummy_array[:,mask] = EHIsig
        oehisig[:] = dummy_array.copy()
        dummy_array = np.ones(original_shape)*const.FILL_VAL
        dummy_array[:,mask] = EHIaccl
        oehiaccl[:] = dummy_array.copy()
    else:
        oehisig[:] = EHIsig
        oehiaccl[:] = EHIaccl
    tempnc.close()
    dailyout.close()
