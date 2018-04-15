#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
ncio.py contains the functions for reading from and writing to NetCDF4 files.
Created on Sat Apr 14 12:56:42 2018

@author: Tammas Loughran
"""
import sys
import datetime as dt
try:
    modulename = 'netCDF4'
    import netCDF4 as nc
    from netCDF4 import MFDataset, MFTime, Dataset
    modulename = 'pandas'
    import pandas as pd
except ImportError:
    print(modulename, " is missing. Please install missing packages.")
    sys.exit(2)
import numpy as np


class DatesOrderError(Exception):
    """Exception to be raised when the calendar start day occurs before the end day"""

    def __init__(self, start, end):
        self.start = start
        self.end = end
        print("Calendar end date appears before the start date")
        print("start: ", self.start)
        print("end: ", self.end)


class Calendar360():
    """Creates a basic calendar object with 360 days per year."""

    def __init__(self,sdate,edate):
        if sdate>edate: raise DatesOrderError(sdate,edate)
        self.year = np.repeat(range(sdate.year,edate.year+1), 360, 0)
        nyears = len(range(sdate.year,edate.year+1))
        self.month = np.tile(np.repeat(range(1,12+1), 30, 0), nyears)
        self.day = np.tile(np.tile(range(1,30+1), 12), nyears)
        if (sdate.day!=1)|(edate.month!=1):
            sdoyi = (sdate.month-1)*30+sdate.day-1
            self.year = self.year[sdoyi:]
            self.month = self.month[sdoyi:]
            self.day = self.day[sdoyi:]
        if (edate.day!=30)|(edate.month!=12):
            edoyi = (12-edate.month)*30+(30-edate.day)
            self.year = self.year[:-edoyi]
            self.month = self.month[:-edoyi]
            self.day = self.day[:-edoyi]

    def __getitem__(self,i):
        """Return a datetime object for the provided index"""
        return dt.datetime(self.year[i], self.month[i], self.day[i])


class TimeData(object):
    """A class to contain all of the time data from an netcdf file."""


def get_time_data(options):
    """This function fetches the time and calendar data from the input netcdf files.

    The argument is the options object, which contains the tmaxfile filename.

    Returns a TimeData object containing a collection of variables.
    """
    timedata = TimeData()
    if any([(wildcard in options.tmaxfile) for wildcard in ['*','?','[']]):
        tmaxnc = MFDataset(options.tmaxfile, 'r')
        nctime = tmaxnc.variables[options.timevname]
        nctime = MFTime(nctime)
    else:
        tmaxnc = Dataset(options.tmaxfile, 'r')
        nctime = tmaxnc.variables[options.timevname]
    timedata.calendar = nctime.calendar

    if not timedata.calendar:
        print('Unrecognized calendar. Using gregorian.')
        timedata.calendar = 'gregorian'

    if timedata.calendar=='360_day':
        timedata.daysinyear = 360
        # 360 day season start and end indices
        timedata.SHS = (300,450)
        timedata.SHW = (120,270)
        timedata.dayone = nc.num2date(nctime[0], nctime.units, calendar=timedata.calendar)
        timedata.daylast = nc.num2date(nctime[-1], nctime.units, calendar=timedata.calendar)
        timedata.dates = Calendar360(timedata.dayone, timedata.daylast)
    else:
        timedata.daysinyear = 365
        # 365 day season start and end indices
        timedata.SHS = (304,455)
        timedata.SHW = (120,273)
        if tmaxnc.variables[options.timevname].units=='day as %Y%m%d.%f':
            st = str(int(nctime[0]))
            nd = str(int(nctime[-1]))
            timedata.dayone = dt.datetime(int(st[:4]), int(st[4:6]), int(st[6:]))
            timedata.daylast = dt.datetime(int(nd[:4]), int(nd[4:6]), int(nd[6:]))
        else:
            timedata.dayone = nc.num2date(nctime[0], nctime.units, calendar=timedata.calendar)
            timedata.daylast = nc.num2date(nctime[-1], nctime.units, calendar=timedata.calendar)
        timedata.dates = pd.period_range(str(timedata.dayone), str(timedata.daylast))
        # Remove leap days. Maybe this should be a separate function?
        if timedata.calendar=='365_day': timedata.dates = timedata.dates[(timedata.dates.month!=2)|(timedata.dates.day!=29)]

        return timedata


def get_mask(options):
    """get_mask loads the land-sea mask data from a NetCDF4 file."""
    masknc = Dataset(options.maskfile, 'r')
    mask = masknc.variables[options.maskvname][:]
    if mask.max()>1: mask = mask>50
    mask = mask.astype(np.bool)
    mask = np.squeeze(mask)
    masknc.close()
    return mask


def load_bp_data(options, timedata, variable='tmax', mask=None):
    """load_bp_data loads the tmax or tmin data for the baseperiod provided in options."""
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
        tempnc = MFDataset(files, 'r')
        bptime = tempnc.variables[options.timevname]
        bptime = MFTime(bptime)
    else:
        tempnc = Dataset(files,'r')
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

    temp = tempnc.variables[varname][(options.bpstart<=bpdates.year)&(bpdates.year<=options.bpend)]
    if len(temp.shape)==4: temp = temp.squeeze()

    # Test for increasing latitude and flip if decreasing
    try:
        lats = tempnc.variables['lat'][:]
    except KeyError:
        lats = tempnc.variables['latitude'][:]
    increasing = (lats[0]-lats[-1])<0
    if not increasing:
        lats = np.flipud(lats)
        temp = np.fliplr(temp)
    if options.maskfile:
        if not increasing: mask = np.flipud(mask)
        temp = temp[:,mask]

    if tempnc.variables[options.tmaxvname].units=='K': temp -= 273.15

    # Remove leap days in gregorian calendars
    if (timedata.calendar=='gregorian')|(timedata.calendar=='proleptic_gregorian')|(timedata.calendar=='standard'):
        temp = temp[(dates_base.month!=2)|(dates_base.day!=29),...]
        del dates_base

    return temp


def save_yearly(HWA,HWM,HWN,HWF,HWD,HWT,tpct,definition):
    """Save yearly data to netcdf file."""
    yearlyout = Dataset('%s_heatwaves_%s_%s_%s_yearly_%s.nc'%(definition,
        model, experiment, rip, season), 'w')
    yearlyout.createDimension('time', size=None)
    yearlyout.createDimension('lon', tmaxnc.dimensions[lonname].__len__())
    yearlyout.createDimension('lat', tmaxnc.dimensions[latname].__len__())
    yearlyout.createDimension('day', daysinyear)
    setattr(yearlyout, "author", "Tammas Loughran")
    setattr(yearlyout, "contact", "t.loughran@student.unsw.edu.au")
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
    setattr(yearlyout, "period", "%s-%s"%(str(first_year),str(daylast.year)))
    setattr(yearlyout, "base_period", "%s-%s"%(str(bpstart),str(bpend)))
    setattr(yearlyout, "percentile", "%sth"%(str(pcntl)))
    setattr(yearlyout, "definition", definition)
    setattr(yearlyout, "frequency", "yearly")
    setattr(yearlyout, "season", season)
    setattr(yearlyout, "season_note", ("The year of a season is the year it starts"
            " in. SH summer: Nov-Mar. NH summer: May-Sep."))
    try:
        file = open('version', 'r')
        commit = file.read()[:]
        if commit[-2:]==r'\n': commit = commit[:-2]
    except IOError:
        commit = "Unknown. Check date for latest version."
    setattr(yearlyout, "git_commit", commit)
    setattr(yearlyout, "tmax_file", options.tmaxfile)
    setattr(yearlyout, "tmin_file", options.tminfile)
    if options.maskfile:
        setattr(yearlyout, "mask_file", options.maskfile)
    setattr(yearlyout, "quantile_method", options.qtilemethod)
    otime = yearlyout.createVariable('time', 'f8', 'time',
            fill_value=-999.99)
    setattr(otime, 'units', 'year')
    olat = yearlyout.createVariable('lat', 'f8', 'lat')
    setattr(olat, 'standard_name', 'latitude')
    setattr(olat, 'long_name', 'Latitude')
    setattr(olat, 'units', 'degrees_north')
    setattr(olat, 'axis', 'Y')
    olon = yearlyout.createVariable('lon', 'f8', 'lon')
    setattr(olon, 'standard_name', 'longiitude')
    setattr(olon, 'long_name', 'Longitude')
    setattr(olon, 'units', 'degrees_east')
    setattr(olon, 'axis', 'X')
    otpct = yearlyout.createVariable('t%spct'%(pcntl), 'f8',
	    ('day','lat','lon'), fill_value=-999.99)
    setattr(otpct, 'long_name', '90th percentile')
    setattr(otpct, 'units', 'degC')
    setattr(otpct, 'description',
            '90th percentile of %s-%s'%(str(bpstart),str(bpend)))
    HWAout = yearlyout.createVariable('HWA_%s'%(definition), 'f8',
            ('time','lat','lon'), fill_value=-999.99)
    setattr(HWAout, 'long_name', 'Heatwave Amplitude')
    if definition=='tx90pct':
        setattr(HWAout, 'units', 'degC')
    elif definition=='tn90pct':
        setattr(HWAout, 'units', 'degC')
    elif definition=='EHF':
        setattr(HWAout, 'units', 'degC2')
    setattr(HWAout, 'description',
            'Peak of the hottest heatwave per year')
    HWMout = yearlyout.createVariable('HWM_%s'%(definition), 'f8',
            ('time','lat','lon'), fill_value=-999.99)
    setattr(HWMout, 'long_name', 'Heatwave Magnitude')
    if definition=='tx90pct':
        setattr(HWAout, 'units', 'degC')
    elif definition=='tn90pct':
        setattr(HWAout, 'units', 'degC')
    elif definition=='EHF':
        setattr(HWAout, 'units', 'degC2')
    setattr(HWMout, 'description', 'Average magnitude of the yearly heatwave')
    HWNout = yearlyout.createVariable('HWN_%s'%(definition), 'f8',
            ('time', 'lat', 'lon'), fill_value=-999.99)
    setattr(HWNout, 'long_name', 'Heatwave Number')
    setattr(HWNout, 'units','heatwaves')
    setattr(HWNout, 'description', 'Number of heatwaves per year')
    HWFout = yearlyout.createVariable('HWF_%s'%(definition), 'f8',
            ('time','lat','lon'), fill_value=-999.99)
    setattr(HWFout, 'long_name','Heatwave Frequency')
    setattr(HWFout, 'units', 'days')
    setattr(HWFout, 'description', 'Proportion of heatwave days per season')
    HWDout = yearlyout.createVariable('HWD_%s'%(definition), 'f8',
            ('time','lat','lon'), fill_value=-999.99)
    setattr(HWDout, 'long_name', 'Heatwave Duration')
    setattr(HWDout, 'units', 'days')
    setattr(HWDout, 'description', 'Duration of the longest heatwave per year')
    HWTout = yearlyout.createVariable('HWT_%s'%(definition), 'f8',
            ('time','lat','lon'), fill_value=-999.99)
    setattr(HWTout, 'long_name', 'Heatwave Timing')
    setattr(HWTout, 'units', 'days from strat of season')
    setattr(HWTout, 'description', 'First heat wave day of the season')
    otime[:] = range(first_year, daylast.year+1)
    olat[:] = lats
    lons = tmaxnc.variables[lonname][:]
    olon[:] = lons
    dummy_array = np.ones((daysinyear,)+(len(lats),)+(len(lons),))*np.nan
    if options.maskfile:
        dummy_array[:,mask] = tpct
        dummy_array[np.isnan(dummy_array)] = -999.99
        otpct[:] = dummy_array.copy()
        dummy_array = np.ones((nyears,)+(len(lats),)+(len(lons),))*np.nan
        dummy_array[:,mask] = HWA
        dummy_array[np.isnan(dummy_array)] = -999.99
        HWAout[:] = dummy_array.copy()
        dummy_array[:,mask] = HWM
        dummy_array[np.isnan(dummy_array)] = -999.99
        HWMout[:] = dummy_array.copy()
        dummy_array[:,mask] = HWN
        dummy_array[np.isnan(dummy_array)] = -999.99
        HWNout[:] = dummy_array.copy()
        dummy_array[:,mask] = HWF
        dummy_array[np.isnan(dummy_array)] = -999.99
        HWFout[:] = dummy_array.copy()
        dummy_array[:,mask] = HWD
        dummy_array[np.isnan(dummy_array)] = -999.99
        HWDout[:] = dummy_array.copy()
        dummy_array[:,mask] = HWT
        dummy_array[np.isnan(dummy_array)] = -999.99
        HWTout[:] = dummy_array.copy()
    else:
        otpct[:] = tpct
        HWAout[:] = HWA.reshape((nyears,)+space)
        HWMout[:] = HWM.reshape((nyears,)+space)
        HWNout[:] = HWN.reshape((nyears,)+space)
        HWFout[:] = HWF.reshape((nyears,)+space)
        HWDout[:] = HWD.reshape((nyears,)+space)
        HWTout[:] = HWT.reshape((nyears,)+space)
    yearlyout.close()