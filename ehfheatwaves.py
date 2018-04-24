#!/usr/bin/python3
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
    modulename = 'numpy'
    import numpy as np
    modulename = 'distutils.version'
    from distutils.version import LooseVersion
except ImportError:
    print(modulename, " is missing. Please install missing packages.")
    sys.exit(2)
if LooseVersion(np.__version__) < LooseVersion('1.8.0'):
    print("Please install numpy version 1.8.0 or higher.")
    sys.exit(2)
import qtiler
import getoptions
import ncio


def window_percentile(temp, options, daysinyear=365, wsize=15):
    """window_percentile calculates a day-of-year moving window percentile."""
    # Initialise array.
    pctl = np.ones(((daysinyear,)+tmax.shape[1:]))*np.nan

    # Construct the window.
    window = np.zeros(daysinyear,dtype=np.bool)
    window[-np.int(np.floor(wsize/2.)):] = 1
    window[:np.int(np.ceil(wsize/2.))] = 1
    window = np.tile(window,options.bpend+1-options.bpstart)

    # Select the interpolation method.
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

    # Set the percentile for each day of year.
    for day in range(daysinyear):
        pctl[day,...] = percentile(temp[window,...], options.pcntl, parameter)
        window = np.roll(window,1)

    return pctl


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


def hw_aspects(EHF, season, hemisphere):
    """hw_aspects takes EHF values or temp 90pct exceedences identifies
    heatwaves and calculates seasonal aspects.
    """
    # Select indices depending on calendar season and hemisphere
    if season=='summer':
        if hemisphere=='south':
            startday = timedata.SHS[0]
            endday = timedata.SHS[1]
        else:
            startday = timedata.SHW[0]
            endday = timedata.SHW[1]
    elif season=='winter':
        if hemisphere=='south':
            startday = timedata.SHW[0]
            endday = timedata.SHW[1]
        else:
            startday = timedata.SHS[0]
            endday = timedata.SHS[1]
    # Initialize arrays
    HWA = np.ones(((nyears,)+(EHF.shape[1],)))*np.nan
    HWM = HWA.copy()
    HWN = HWA.copy()
    HWF = HWA.copy()
    HWD = HWA.copy()
    HWT = HWA.copy()
    # Loop over years
    for iyear, year in enumerate(range(first_year,timedata.daylast.year)):
        if options.oldmethod:
            if (year==timedata.daylast.year): continue # Incomplete yr
            # Select this years season
            allowance = 14 # For including heawave days after the end of the season
            ifrom = timedata.startday + timedata.daysinyear*iyear - 1 # -1 to include Oct 31st
            ito = endday + timedata.daysinyear*iyear + allowance
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
            ifrom = startday + timedata.daysinyear*iyear - 2
            ito = endday + timedata.daysinyear*iyear + allowance
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

    # Load time data
    if options.verbose: print("Loading data")
    timedata = ncio.get_time_data(options)

    # Load land-sea mask
    if options.maskfile: mask = ncio.get_mask(options)
    else: mask = None

    # Load the temperature data over the base period
    if options.keeptave or options.keeptmax:
        tmax = ncio.load_bp_data(options, timedata, variable='tmax', mask=mask)
    if options.keeptave or options.keeptmin:
        tmin = ncio.load_bp_data(options, timedata, variable='tmin', mask=mask)
    if options.keeptave:
        tave_base = (tmax + tmin)/2.

    # Caclulate percentile
    if options.verbose: print("Calculating percentiles")
    if options.keeptave:
        tpct = window_percentile(tave_base, options, daysinyear=timedata.daysinyear)
    if options.keeptmax:
        txpct = window_percentile(tmax, options, daysinyear=timedata.daysinyear)
    if options.keeptmin:
        tnpct = window_percentile(tmin, options, daysinyear=timedata.daysinyear)
    if not options.noehf: del tave_base

    # Load all data
    if options.verbose:print("Loading data")
    if options.keeptave or options.keeptmax:
        tmax, lats = ncio.get_all_data(options.tmaxfile, options.tmaxvname, options)
        original_shape = tmax.shape
    if options.keeptave or options.keeptmin:
        tmin, _ = ncio.get_all_data(options.tminfile, options.tminvname, options)
        original_shape = tmin.shape
    if options.keeptave:
        tave = (tmax + tmin)/2.

    # Remove leap days from tave
    if (timedata.calendar=='gregorian')|(timedata.calendar=='proleptic_gregorian')|(timedata.calendar=='standard'):
        if options.keeptave:
            tave = tave[(timedata.dates.month!=2)|(timedata.dates.day!=29),...]
            original_shape = (tave.shape[0], original_shape[1], original_shape[2])
        if options.keeptmax or options.keeptave:
            tmax = tmax[(timedata.dates.month!=2)|(timedata.dates.day!=29),...]
            original_shape = (tmax.shape[0], original_shape[1], original_shape[2])
        if options.keeptmin or options.keeptave:
            tmin = tmin[(timedata.dates.month!=2)|(timedata.dates.day!=29),...]
            original_shape = (tmin.shape[0], original_shape[1], original_shape[2])
        calendar = '365_day'

    # Remove incomplete starting year
    first_year = timedata.dayone.year
    if (timedata.dayone.month!=1)|(timedata.dayone.day!=1):
        first_year = timedata.dayone.year+1
        start = np.argmax(timedata.dates.year==first_year)
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
                    tpct[i-timedata.daysinyear*int((i+1)/timedata.daysinyear),...]
            EHF[i,...] = np.maximum(EHIaccl,1.)*EHIsig
        EHF[EHF<0] = 0

    # Tx90pc exceedences
    if options.keeptmin or options.keeptmax:
        if options.keeptmax:
            txexceed = np.ones(tmax.shape)*np.nan
            for i in range(0,tmax.shape[0]):
                idoy = i-timedata.daysinyear*int((i+1)/timedata.daysinyear)
                txexceed[i,...] = tmax[i,...]>txpct[idoy,...]
            txexceed[txexceed>0] = tmax[txexceed>0]
        if options.keeptmin:
            tnexceed = np.ones(tmin.shape)*np.nan
            for i in range(0,tmin.shape[0]):
                idoy = i-timedata.daysinyear*int((i+1)/timedata.daysinyear)
                tnexceed[i,...] = tmin[i,...]>tnpct[idoy,...]
            tnexceed[tnexceed>0] = tmin[tnexceed>0]

    # Calculate daily output
    if options.daily: event, ends = identify_hw(EHF)
    if options.tx90pcd:
        event_tx, ends_tx = identify_hw(txexceed)
    if options.tn90pcd:
        event_tn, ends_tn = identify_hw(tnexceed)

    nyears = len(range(first_year,timedata.daylast.year+1))

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
    # Save yearly data to netcdf
    if options.yearlyout:
        if not options.noehf:
            ncio.save_yearly(HWA_EHF,HWM_EHF,HWN_EHF,HWF_EHF,HWD_EHF,HWT_EHF,tpct,"EHF",timedata,options,mask)
        if options.tx90pc:
            ncio.save_yearly(HWA_tx,HWM_tx,HWN_tx,HWF_tx,HWD_tx,HWT_tx,txpct,"tx90pct",timedata,options,mask)
        if options.tn90pc:
            ncio.save_yearly(HWA_tn,HWM_tn,HWN_tn,HWF_tn,HWD_tn,HWT_tn,tnpct,"tn90pct",timedata,options,mask)

    # Save daily data to netcdf
    if options.dailyout:
        if options.keeptave:
            ncio.save_daily(EHF, event, ends, options, options, timedata, original_shape, mask=mask)
        if options.tx90pcd:
            ncio.save_daily(txexceed, event_tx, ends_tx, options, timedata, original_shape, mask=mask)
        if options.tn90pcd:
            ncio.save_daily(tnexceed, event_tn, ends_tn, options, timedata, original_shape, mask=mask)