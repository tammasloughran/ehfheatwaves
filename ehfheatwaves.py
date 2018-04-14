#!/usr/bin/ python3
# -*- coding: utf-8 -*-
"""
ehfheatwaves.py calculates heatwave indices and characteristics from
temperature data.

@author: Tammas Loughran
"""

import sys
import warnings
warnings.simplefilter('ignore',category=RuntimeWarning)
try:
    modulename = 'pandas'
    import pandas as pd
    modulename = 'numpy'
    import numpy as np
    modulename = 'datetime'
    import datetime as dt
    modulename = 'qtiler'
    import qtiler
    modulename = 'netCDF4'
    import netCDF4 as nc
    from netCDF4 import MFDataset, MFTime, Dataset
    modulename = 'distutils.version'
    from distutils.version import LooseVersion
except ImportError:
    print(modulename, " is missing. Please install missing packages.")
    sys.exit(2)
if LooseVersion(np.__version__) < LooseVersion('1.8.0'):
    print("Please install numpy version 1.8.0 or higher.")
    sys.exit(2)
import getoptions


class DatesOrderError(Exception):
    """Exception to be raised when the calendar start day occurs before the end day"""

    def __init__(self, start, end):
        self.start = start
        self.end = end
        print("Calendar end date appears before the start date")
        print("start: ", self.start)
        print("end: ", self.end)


def identify_hw(ehfs):
    """identify_hw locates heatwaves from EHF and returns an event indicator
    and a duration indicator.
    """
    # Agregate consecutive days with EHF>0
    # First day contains duration
    events = (ehfs!=0.).astype(np.int)
    for i in range(events.shape[0]-2,-1,-1):
         events[i,events[i,...]>0] = events[i+1,events[i,...]>0]+1

    # Identify when heatwaves start with duration
    # Given that first day contains duration
    diff = np.zeros(events.shape)
    # Insert the first diff value as np.diff doesn't catch it because
    # there is no pevious value to compare to.
    diff[0,...] = events[0,...]
    diff[1:,...] = np.diff(events, axis=0)
    endss = np.zeros(ehfs.shape,dtype=np.int)
    endss[diff>2] = events[diff>2]

    # Remove events less than 3 days
    events[diff==2] = 0
    events[np.roll(diff==2, 1, axis=0)] = 0
    events[diff==1] = 0
    del diff
    events[events>0] = 1
    events = events.astype(np.bool)
    endss[endss<3] = 0
    return events, endss


def identify_semi_hw(ehfs):
    """identify_hw locates heatwaves from EHF and returns an event indicator
    and a duration indicator. This function does not exclude events less than
    three days in duration.
    """
    # Agregate consecutive days with EHF>0
    # First day contains duration
    events = (ehfs!=0.).astype(np.int)
    for i in range(events.shape[0]-2,-1,-1):
         events[i,events[i,...]>0] = events[i+1,events[i,...]>0]+1

    # Identify when heatwaves start with duration
    # Given that first day contains duration
    diff = np.zeros(events.shape)
    # Insert the first diff value as np.diff doesn't catch it because
    # there is no pevious value to compare to.
    diff[0,...] = events[0,...]
    diff[1:,...] = np.diff(events, axis=0)
    endss = np.zeros(ehfs.shape,dtype=np.int)
    endss[diff>0] = events[diff>0]

    del diff
    events[events>0] = 1
    events = events.astype(np.bool)
    return events, endss


class calendar360():
    """Creates a calendar object with 360 days per year."""

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


def hw_aspects(EHF, season, hemisphere):
    """hw_aspects takes EHF values or temp 90pct exceedences identifies
    heatwaves and calculates seasonal aspects.
    """
    # Select indices depending on calendar season and hemisphere
    if season=='summer':
        if hemisphere=='south':
            startday = SHS[0]
            endday = SHS[1]
        else:
            startday = SHW[0]
            endday = SHW[1]
    elif season=='winter':
        if hemisphere=='south':
            startday = SHW[0]
            endday = SHW[1]
        else:
            startday = SHS[0]
            endday = SHS[1]
    # Initialize arrays
    HWA = np.ones(((nyears,)+(EHF.shape[1],)))*np.nan
    HWM = HWA.copy()
    HWN = HWA.copy()
    HWF = HWA.copy()
    HWD = HWA.copy()
    HWT = HWA.copy()
    # Loop over years
    for iyear, year in enumerate(range(first_year,daylast.year)):
        if options.oldmethod:
            if (year==daylast.year): continue # Incomplete yr
            # Select this years season
            allowance = 14 # For including heawave days after the end of the season
            ifrom = startday + daysinyear*iyear - 1 # -1 to include Oct 31st
            ito = endday + daysinyear*iyear + allowance
            EHF_i = EHF[ifrom:ito,...]
            event_i, duration_i = identify_hw(EHF_i)
            # Identify heatwaves that span the entire season
            perpetual = event_i[:-allowance,...].all(axis=0)
            perphw = duration_i[0,perpetual] - 1 # -1 to exclude Oct 31st
            # Remove events that start after the end of the season and before start
            EHF_i = EHF_i[1:,...]
            duration_i = duration_i[1:-allowance,...]
            event_i = event_i[1:-allowance,...]
            # Indicate perpetual heatwaves if they occur.
            if perpetual.any(): duration_i[0,perpetual] = perphw
        else:
            # Select this years season
            allowance = 14 # For including heawave days after the end of the season
            ifrom = startday + daysinyear*iyear - 2
            ito = endday + daysinyear*iyear + allowance
            EHF_i = EHF[ifrom:ito,...]
            event_i, duration_i = identify_hw(EHF_i)
            # Remove EHF values in pre season
            EHF_i = EHF_i[2:,...]
            # Identify semi heatwaves that overlap the start of season only including days within the season.
            event_i = event_i[2:,...]
            event_i, duration_i = identify_semi_hw(event_i)
            # Identify heatwaves that span the entire season
            perpetual = event_i[:-allowance,...].all(axis=0)
            perphw = duration_i[0,perpetual] - 1 # -1 to exclude Oct 31st
            # Indicate locations of perpetual heatwaves if they occur.
            if perpetual.any(): duration_i[0,perpetual] = perphw
            # Remove events that start after the end of the season
            duration_i = duration_i[:-allowance,...]
        # Calculate metrics
        HWN[iyear,...] = (duration_i>0).sum(axis=0)
        HWF[iyear,...] = duration_i.sum(axis=0)
        HWD[iyear,...] = duration_i.max(axis=0)
        HWT[iyear,...] = np.argmax(event_i,axis=0)
        HWT[iyear,HWD[iyear,...]==0] = np.nan
        HWD[iyear,HWD[iyear,...]==0] = np.nan
        # HWM and HWA must be done on each gridcell
        for x in range(EHF_i.shape[1]):
            hw_mag = []
            # retrieve indices where heatwaves start.
            i = np.where(duration_i[:,x]>0)[0] # time
            d = duration_i[i,x] # duration
            if (d==0).all(): continue
            for hw in range(len(d)):
                # retireve this heatwave's EHF values and mean magnitude
                hwdat = EHF_i[i[hw]:i[hw]+d[hw],x]
                hw_mag.append(np.nanmean(hwdat))
            HWM[iyear,x] = np.nanmean(hw_mag)
            # Find the hottest heatwave magnitude
            idex = np.where(hw_mag==max(hw_mag))[0][0]
            # Find that heatwave's hottest day as EHF value.
            HWA[iyear,x] = EHF_i[i[idex]:i[idex]+d[idex],x].max()
    return HWA, HWM, HWN, HWF, HWD, HWT


# Calculate metrics year by year
def split_hemispheres(EHF):
    """split_hemispheres splits the input data by hemispheres, and glues them
    back together after heatwave calculations.
    """
    if south:
        if options.maskfile:
            EHF_s = EHF[:,:(mask[lats<=0]>0).sum()]
        else:
            EHF_s = EHF[:,lats<=0,...]
        # Reshape to 2D
        space = EHF_s.shape[1:]
        if len(space)>1:
            EHF_s = EHF_s.reshape(EHF_s.shape[0],space[0]*space[1])
        # Southern hemisphere aspects
        HWA_s, HWM_s, HWN_s, HWF_s, HWD_s, HWT_s = \
                hw_aspects(EHF_s, options.season, 'south')
        del EHF_s
    if north:
        if options.maskfile:
            EHF_n = EHF[:,(mask[lats<=0]>0).sum():]
        else:
            EHF_n = EHF[:,lats>0,...]
        # Reshape to 2D
        space = EHF_n.shape[1:]
        if len(space)>1:
            EHF_n = EHF_n.reshape(EHF_n.shape[0],space[0]*space[1])
        # Northern hemisphere aspects
        HWA_n, HWM_n, HWN_n, HWF_n, HWD_n, HWT_n = \
                hw_aspects(EHF_n, options.season, 'north')
        del EHF_n
    # Glue hemispheres back together
    if north and south:
        HWA = np.append(HWA_s, HWA_n, axis=1)
        HWM = np.append(HWM_s, HWM_n, axis=1)
        HWN = np.append(HWN_s, HWN_n, axis=1)
        HWF = np.append(HWF_s, HWF_n, axis=1)
        HWD = np.append(HWD_s, HWD_n, axis=1)
        HWT = np.append(HWT_s, HWT_n, axis=1)
    elif north:
        HWA = HWA_n
        HWM = HWM_n
        HWN = HWN_n
        HWF = HWF_n
        HWD = HWD_n
        HWT = HWT_n
    elif south:
        HWA = HWA_s
        HWM = HWM_s
        HWN = HWN_s
        HWF = HWF_s
        HWD = HWD_s
        HWT = HWT_s
    return HWA, HWM, HWN, HWF, HWD, HWT

if __name__=='__main__':

    # Get the options and variables
    options = getoptions.parse_arguments(sys.argv[1:])

    if options.verbose: print("Loading data")
    # Load time data
    try:
        tmaxnc = MFDataset(options.tmaxfile, 'r')
    except IndexError:
        tmaxnc = Dataset(options.tmaxfile, 'r')
    nctime = tmaxnc.variables[options.timevname]
    try:
        nctime = MFTime(nctime)
    except AttributeError:
        pass
    except ValueError:
        pass
    calendar = nctime.calendar
    if not calendar:
        print('Unrecognized calendar. Using gregorian.')
        calendar = 'gregorian'
    elif calendar=='360_day':
        daysinyear = 360
        # 360 day season start and end indices
        SHS = (300,450)
        SHW = (120,270)
        dayone = nc.num2date(nctime[0], nctime.units,
                calendar=calendar)
        daylast = nc.num2date(nctime[-1], nctime.units,
                calendar=calendar)
        dates = calendar360(dayone, daylast)
        shorten = 0
        if (daylast.day!=30)|(daylast.month!=12):
            shorten = 30*(13-daylast.month) - daylast.day
    else:
        daysinyear = 365
        # 365 day season start and end indices
        SHS = (304,455)
        SHW = (120,273)
        if tmaxnc.variables[options.timevname].units=='day as %Y%m%d.%f':
            st = str(int(nctime[0]))
            nd = str(int(nctime[-1]))
            dayone = dt.datetime(int(st[:4]), int(st[4:6]), int(st[6:]))
            daylast = dt.datetime(int(nd[:4]), int(nd[4:6]), int(nd[6:]))
        else:
            dayone = nc.num2date(nctime[0], nctime.units,
                    calendar=calendar)
            daylast = nc.num2date(nctime[-1], nctime.units,
                    calendar=calendar)
        dates = pd.period_range(str(dayone), str(daylast))
        if calendar=='365_day': dates = dates[(dates.month!=2)|(dates.day!=29)]
        shorten = 0
        if (daylast.day!=30)|(daylast.month!=12):
            endofdata = dt.datetime(2000, daylast.month, daylast.day)
            shorten = dt.datetime(2000, 12, 31) - endofdata
            shorten = shorten.days

    # Load land-sea mask
    if options.maskfile:
        masknc = Dataset(options.maskfile, 'r')
        vname = options.maskvname
        mask = masknc.variables[vname][:]
        if mask.max()>1: mask = mask>50
        mask = mask.astype(np.bool)
        mask = np.squeeze(mask)
        masknc.close()

    # Load base period data
    if options.bpfn:
        try:
            tminnc = MFDataset(options.bpfn, 'r')
        except IndexError:
            tminnc = Dataset(options.bpfn, 'r')
    else:
        try:
            tminnc = MFDataset(options.tminfile, 'r')
        except IndexError:
            tminnc = Dataset(options.tminfile, 'r')
    if options.bpfx:
        try:
            tmaxnc = MFDataset(options.bpfx, 'r')
        except IndexError:
            tmaxnc = Dataset(options.bpfx, 'r')
    else:
        try:
            tmaxnc = MFDataset(options.tmaxfile, 'r')
        except IndexError:
            tmaxnc = Dataset(options.tmaxfile, 'r')
    vname = options.tmaxvname
    bptime = tmaxnc.variables[options.timevname]
    try:
        bptime = MFTime(bptime)
    except AttributeError:
        pass
    except ValueError:
        pass
    if tmaxnc.variables[options.timevname].units=='day as %Y%m%d.%f':
        st = str(int(bptime[0]))
        nd = str(int(bptime[-1]))
        bpdayone = dt.datetime(int(st[:4]), int(st[4:6]), int(st[6:]))
        bpdaylast = dt.datetime(int(nd[:4]), int(nd[4:6]), int(nd[6:]))
    else:
        bpdayone = nc.num2date(bptime[0], bptime.units, calendar=calendar)
        bpdaylast = nc.num2date(bptime[-1], bptime.units, calendar=calendar)
    if calendar=='360_day': bpdates = calendar360(bpdayone, bpdaylast)
    else:
        bpdates = pd.period_range(str(bpdayone), str(bpdaylast))
        if calendar=='365_day': bpdates = bpdates[(bpdates.month!=2)|(bpdates.day!=29)]
        dates_base = bpdates[(options.bpstart<=bpdates.year)&(bpdates.year<=options.bpend)]
    tmax = tmaxnc.variables[options.tmaxvname][(options.bpstart<=bpdates.year)&(bpdates.year<=options.bpend)]
    tmin = tminnc.variables[options.tminvname][(options.bpstart<=bpdates.year)&(bpdates.year<=options.bpend)]
    if len(tmin.shape)==4: tmin = tmin.squeeze()
    if len(tmax.shape)==4: tmax = tmax.squeeze()
    # Test for increasing latitude and flip if decreasing
    try:
        lats = tmaxnc.variables['lat'][:]
    except KeyError:
        lats = tmaxnc.variables['latitude'][:]
    increasing = (lats[0]-lats[-1])<0
    if not increasing:
        lats = np.flipud(lats)
        tmax = np.fliplr(tmax)
        tmin = np.fliplr(tmin)
    if options.maskfile:
        if not increasing: mask = np.flipud(mask)
        tmin = tmin[:,mask]
        tmax = tmax[:,mask]
    if tmaxnc.variables[options.tmaxvname].units=='K':
        tmax -= 273.15
        tmin -= 273.15
    if options.keeptave: tave_base = (tmax + tmin)/2.
    if not options.keeptmin: del tmin
    if not options.keeptmax: del tmax

    # Remove leap days in gregorian calendars
    if (calendar=='gregorian')|(calendar=='proleptic_gregorian')|\
                (calendar=='standard'):
        if options.keeptave:
            tave_base = tave_base[(dates_base.month!=2)|(dates_base.day!=29),...]
        if options.keeptmax:
            tmax = tmax[(dates_base.month!=2)|(dates_base.day!=29),...]
        if options.keeptmin:
            tmin = tmin[(dates_base.month!=2)|(dates_base.day!=29),...]
        del dates_base

    if options.verbose: print("Calculating percentiles")
    # Caclulate 90th percentile
    if not options.noehf: tpct = np.ones(((daysinyear,)+tave_base.shape[1:]))*np.nan
    if options.keeptmax: txpct = np.ones(((daysinyear,)+tmax.shape[1:]))*np.nan
    if options.keeptmin: tnpct = np.ones(((daysinyear,)+tmin.shape[1:]))*np.nan
    window = np.zeros(daysinyear,dtype=np.bool)
    wsize = 15.
    window[-np.int(np.floor(wsize/2.)):] = 1
    window[:np.int(np.ceil(wsize/2.))] = 1
    window = np.tile(window,options.bpend+1-options.bpstart)
    if options.qtilemethod=='python':
        percentile = np.percentile
        parameter = 0
    elif options.qtilemethod=='zhang':
        percentile = qtiler.quantile_zhang
        parameter = False
    elif options.qtilemethod=='matlab':
        percentile = qtiler.quantile_R
        parameter = 5
    elif options.qtilemethod=='climpact':
        percentile = qtiler.quantile_climpact
        parameter = False
    for day in range(daysinyear):
        if options.keeptave:
            tpct[day,...] = percentile(tave_base[window,...], options.pcntl, parameter)
        if options.keeptmax:
            txpct[day,...] = percentile(tmax[window,...], options.pcntl, parameter)
        if options.keeptmin:
            tnpct[day,...] = percentile(tmin[window,...], options.pcntl, parameter)
        window = np.roll(window,1)
    if not options.noehf: del tave_base
    del window

    if options.verbose:print("Loading data")
    # Load all data
    try:
        tminnc = MFDataset(options.tminfile, 'r')
    except IndexError:
        tminnc = Dataset(options.tminfile, 'r')
    try:
        tmaxnc = MFDataset(options.tmaxfile, 'r')
    except IndexError:
        tmaxnc = Dataset(options.tmaxfile, 'r')
    tmax = tmaxnc.variables[options.tmaxvname][:]
    if len(tmax.shape)==4: tmax = tmax.squeeze()
    original_shape = tmax.shape
    tmin = tminnc.variables[options.tminvname][:]
    if len(tmin.shape)==4: tmin = tmin.squeeze()

    # Test for increasing latitude and flip if decreasing
    try:
        lats = tmaxnc.variables['lat'][:]
    except KeyError:
        lats = tmaxnc.variables['latitude'][:]
    increasing = (lats[0]-lats[-1])<0
    if not increasing:
        lats = np.flipud(lats)
        tmax = np.fliplr(tmax)
        tmin = np.fliplr(tmin)
    if options.maskfile:
    #    if not increasing: mask = np.flipud(mask)
        tmax = tmax[:,mask]
        tmin = tmin[:,mask]
    if tmaxnc.variables[options.tmaxvname].units=='K': tmax -= 273.15
    if tminnc.variables[options.tminvname].units=='K': tmin -= 273.15
    if options.keeptave: tave = (tmax + tmin)/2.
    if not options.keeptmin: del tmin
    if not options.keeptmax: del tmax

    # Remove leap days from tave
    if (calendar=='gregorian')|(calendar=='proleptic_gregorian')|\
                (calendar=='standard'):
        if options.keeptave:
            tave = tave[(dates.month!=2)|(dates.day!=29),...]
            original_shape = (tave.shape[0], original_shape[1], original_shape[2])
        if options.keeptmax:
            tmax = tmax[(dates.month!=2)|(dates.day!=29),...]
            original_shape = (tmax.shape[0], original_shape[1], original_shape[2])
        if options.keeptmin:
            tmin = tmin[(dates.month!=2)|(dates.day!=29),...]
            original_shape = (tmin.shape[0], original_shape[1], original_shape[2])
        calendar = '365_day'

    # Remove incomplete starting year
    first_year = dayone.year
    if (dayone.month!=1)|(dayone.day!=1):
        first_year = dayone.year+1
        start = np.argmax(dates.year==first_year)
        if options.keeptave:
            tave = tave[start:,...]
            original_shape = (tave.shape[0], original_shape[1], original_shape[2])
        if options.keeptmax:
            tmax = tmax[start:,...]
            original_shape = (tmax.shape[0], original_shape[1], original_shape[2])
        if options.keeptmin:
            tmin = tmin[start:,...]
            original_shape = (tmin.shape[0], original_shape[1], original_shape[2])

    if options.verbose: print("Caclulating definition")
    # Calculate EHF
    if not options.noehf:
        EHF = np.ones(tave.shape)*np.nan
        for i in range(32,tave.shape[0]):
            EHIaccl = tave[i-2:i+1,...].sum(axis=0)/3. - \
                    tave[i-32:i-2,...].sum(axis=0)/30.
            EHIsig = tave[i-2:i+1,...].sum(axis=0)/3. - \
                    tpct[i-daysinyear*int((i+1)/daysinyear),...]
            EHF[i,...] = np.maximum(EHIaccl,1.)*EHIsig
        EHF[EHF<0] = 0

    # Tx90pc exceedences
    if options.keeptmin or options.keeptmax:
        if options.keeptmax:
            txexceed = np.ones(tmax.shape)*np.nan
            for i in range(0,tmax.shape[0]):
                idoy = i-daysinyear*int((i+1)/daysinyear)
                txexceed[i,...] = tmax[i,...]>txpct[idoy,...]
            txexceed[txexceed>0] = tmax[txexceed>0]
        if options.keeptmin:
            tnexceed = np.ones(tmin.shape)*np.nan
            for i in range(0,tmin.shape[0]):
                idoy = i-daysinyear*int((i+1)/daysinyear)
                tnexceed[i,...] = tmin[i,...]>tnpct[idoy,...]
            tnexceed[tnexceed>0] = tmin[tnexceed>0]

    # Calculate daily output
    if options.daily: event, ends = identify_hw(EHF)
    if options.tx90pcd:
        event_tx, ends_tx = identify_hw(txexceed)
    if options.tn90pcd:
        event_tn, ends_tn = identify_hw(tnexceed)

    nyears = len(range(first_year,daylast.year+1))

    # Calculate yearly output
    if options.yearlyout:
        if options.verbose: print("Calculating yearly aspects")
        # Split by latitude
        north = (lats>0).any()
        south = (lats<=0).any()
        if not options.noehf:
            HWA_EHF, HWM_EHF, HWN_EHF, HWF_EHF, HWD_EHF, HWT_EHF = \
                    split_hemispheres(EHF)
        if options.tx90pc:
            HWA_tx, HWM_tx, HWN_tx, HWF_tx, HWD_tx, HWT_tx = \
                    split_hemispheres(txexceed)
        if options.tn90pc:
            HWA_tn, HWM_tn, HWN_tn, HWF_tn, HWD_tn, HWT_tn = \
                    split_hemispheres(tnexceed)

    if options.verbose: print("Saving")
    # Save to netCDF
    # Retrieve metadata from file
    try:
        experiment = tmaxnc.__getattribute__('experiment')
        model = tmaxnc.__getattribute__('model_id')
        parent = tmaxnc.__getattribute__('parent_experiment_rip')
        realization = tmaxnc.__getattribute__('realization')
        initialization = tmaxnc.__getattribute__('initialization_method')
        physics = tmaxnc.__getattribute__('physics_version')
        rip = 'r'+str(realization)+'i'+str(initialization)+'p'+str(physics)
    except AttributeError:
        experiment = ''
        model = ''
        realization = ''
        rip = ''
        initialization = ''
        physics = ''
    try:
        space = (tmaxnc.dimensions['lat'].__len__(),
                tmaxnc.dimensions['lon'].__len__())
        lonname = 'lon'
        latname = 'lat'
    except KeyError:
        lonname = 'longitude'
        latname = 'latitude'
        space = (tmaxnc.dimensions['latitude'].__len__(),
                tmaxnc.dimensions['longitude'].__len__())


    def save_yearly(HWA,HWM,HWN,HWF,HWD,HWT,tpct,definition):
        """Save yearly data to netcdf file."""
        yearlyout = Dataset('%s_heatwaves_%s_%s_%s_yearly_%s.nc'%(definition,
            model, experiment, rip, options.season), 'w')
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
        setattr(yearlyout, "base_period", "%s-%s"%(str(options.bpstart),str(options.bpend)))
        setattr(yearlyout, "percentile", "%sth"%(str(options.pcntl)))
        setattr(yearlyout, "definition", definition)
        setattr(yearlyout, "frequency", "yearly")
        setattr(yearlyout, "season", options.season)
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
        otpct = yearlyout.createVariable('t%spct'%(options.pcntl), 'f8',
    	    ('day','lat','lon'), fill_value=-999.99)
        setattr(otpct, 'long_name', '90th percentile')
        setattr(otpct, 'units', 'degC')
        setattr(otpct, 'description',
                '90th percentile of %s-%s'%(str(options.bpstart),str(options.bpend)))
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


    # Save yearly data to netcdf
    if options.yearlyout:
        if not options.noehf:
            save_yearly(HWA_EHF,HWM_EHF,HWN_EHF,HWF_EHF,HWD_EHF,HWT_EHF,tpct,"EHF")
        if options.tx90pc:
            save_yearly(HWA_tx,HWM_tx,HWN_tx,HWF_tx,HWD_tx,HWT_tx,txpct,"tx90pct")
        if options.tn90pc:
            save_yearly(HWA_tn,HWM_tn,HWN_tn,HWF_tn,HWD_tn,HWT_tn,tnpct,"tn90pct")


    # Save daily data to netcdf
    if options.dailyout:
        if options.daily or options.dailyonly: defn ='EHF'
        elif options.tx90pcd: defn = 'tx90pct'
        elif options.tn90pcd: defn = 'tn90pct'
        dailyout = Dataset('%s_heatwaves_%s_%s_%s_daily.nc'\
                %(defn, model, experiment, rip), mode='w')
        dailyout.createDimension('time', size=None)
        dailyout.createDimension('lon', tmaxnc.dimensions[lonname].__len__())
        dailyout.createDimension('lat', tmaxnc.dimensions[latname].__len__())
        setattr(dailyout, "author", "Tammas Loughran")
        setattr(dailyout, "contact", "t.loughran@student.unsw.edu.au")
        setattr(dailyout, "source", "https://github.com/tammasloughran/ehfheatwaves")
        setattr(dailyout, "date", dt.datetime.today().strftime('%Y-%m-%d'))
        setattr(dailyout, "script", sys.argv[0])
        setattr(dailyout, "period", "%s-%s"%(str(first_year),str(daylast.year)))
        setattr(dailyout, "base_period", "%s-%s"%(str(options.bpstart),str(options.bpend)))
        setattr(dailyout, "percentile", "%sth"%(str(options.pcntl)))
        if model:
            setattr(dailyout, "model_id", model)
            setattr(dailyout, "experiment", experiment)
            setattr(dailyout, "parent_experiment_rip", parent)
            setattr(dailyout, "realization", realization)
            setattr(dailyout, "initialization_method", initialization)
            setattr(dailyout, "physics_version", physics)
        try:
            file = open('version', 'r')
            commit = file.read()[:]
            if commit[-2:]==r'\n': commit = commit[:-2]
        except IOError:
            commit = "Unknown. Check date for latest version."
        setattr(dailyout, "git_commit", commit)
        setattr(dailyout, "tmax_file", options.tmaxfile)
        setattr(dailyout, "tmin_file", options.tminfile)
        if options.maskfile:
            setattr(dailyout, "mask_file", str(options.maskfile))
        setattr(dailyout, "quantile_method", options.qtilemethod)
        otime = dailyout.createVariable('time', 'f8', 'time',
                        fill_value=-999.99)
        setattr(otime, 'units', 'days since %s-01-01'%(first_year))
        setattr(otime, 'calendar', calendar)
        olat = dailyout.createVariable('lat', 'f8', 'lat')
        setattr(olat, 'standard_name', 'latitude')
        setattr(olat, 'long_name', 'Latitude')
        setattr(olat, 'units', 'degrees_north')
        olon = dailyout.createVariable('lon', 'f8', 'lon')
        setattr(olon, 'standard_name', 'longitude')
        setattr(olon, 'long_name', 'Longitude')
        setattr(olon, 'units', 'degrees_east')
        oehf = dailyout.createVariable(defn, 'f8',
                ('time','lat','lon'), fill_value=-999.99)
        setattr(oehf, 'standard_name', defn)
        if defn=='EHF':
            setattr(oehf, 'long_name', 'Excess Heat Factor')
            setattr(oehf, 'units', 'degC2')
        elif defn=='tx90pct':
            setattr(oehf, 'long_name', 'Temperature Exceeding tx90pct')
            setattr(oehf, 'units', 'C')
        elif defn=='tn90pct':
            setattr(oehf, 'long_name', 'Temperature Exceeding tn90pct')
            setattr(oehf, 'units', 'C')
        oevent = dailyout.createVariable('event', 'f8',
                ('time','lat','lon'), fill_value=-999.99)
        setattr(oevent, 'long_name', 'Event indicator')
        setattr(oevent, 'description',
                'Indicates whether a heatwave is happening on that day')
        oends = dailyout.createVariable('ends', 'f8',
                ('time','lat','lon'), fill_value=-999.99)
        setattr(oends, 'long_name', 'Duration at start of heatwave')
        setattr(oends, 'units', 'days')
        otime[:] = range(0,original_shape[0],1)
        olat[:] = lats
        olon[:] = tmaxnc.variables[lonname][:]
        if options.maskfile:
            dummy_array = np.ones(original_shape)*np.nan
            if options.daily: dummy_array[:,mask] = EHF
            elif options.tx90pcd: dummy_array[:,mask] = txexceed
            elif options.tn90pcd: dummy_array[:,mask] = tnexceed
            dummy_array[np.isnan(dummy_array)] = -999.99
            oehf[:] = dummy_array.copy()
            if options.daily: dummy_array[:,mask] = event
            elif options.tx90pcd: dummy_array[:,mask] = event_tx
            elif options.tn90pcd: dummy_array[:,mask] = event_tn
            dummy_array[np.isnan(dummy_array)] = -999.99
            dummy_array[:31,...] = -999.99
            oevent[:] = dummy_array.copy()
            if options.daily: dummy_array[:,mask] = ends
            elif options.tx90pcd: dummy_array[:,mask] = ends_tx
            elif options.tn90pcd: dummy_array[:,mask] = ends_tn
            dummy_array[np.isnan(dummy_array)] = -999.99
            dummy_array[:31,...] = -999.99
            oends[:] = dummy_array.copy()
        else:
            if options.daily:
                oehf[:] = EHF
                oevent[:] = event
                oends[:] = ends
            elif options.tx90pcd:
                oehf[:] = txexceed
                oevent[:] = event_tx
                oends[:] = ends_tx
            elif options.tn90pcd:
                oehf[:] = tnexceed
                oevent[:] = event_tn
                oends[:] = ends_tn
        dailyout.close()