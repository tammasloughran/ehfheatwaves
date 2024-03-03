# -*- coding: utf-8 -*-
"""Main heatwave calculation module. The main program routine is executed from this module.

The procedure is:

- parse arguments
- load the base period data
- calculate the percentile thresholds
- load the remaining data
- calculate heatwave indices
- calculate seasonal characteristics for each hemisphere and join them back together
- save to netCDF files.

The EHF Index is defined from Nairn et al., (2009) and Nairn and Fawcett, (2013).

$$ EHI_{sig} = \\frac{(T_i + T_{i-1} + T_{i-2})}{3} - T_{90} $$

$$ EHI_{accl} = \\frac{(T_i + T_{i-1} + T_{i-2})}{3} - \\frac{(T_i + ... + T_{-30})}{30} $$

$$ EHF = EHI_{sig} \cdot max(1, EHI_{accl}) $$

Nairn J, Fawcett R. 2013. Defining heatwaves: heatwave defined as a heat-impact event servicing
all community and business sectors in Australia. CAWCR Tech. Rep. 60: 10–15.

Nairn J, Fawcett R, Ray D. 2009. Defining and predicting excessive heat events, a national system.
CAWCR Tech. Rep. 60: 83–86.
"""
import sys
import warnings

import numpy as np

import ehfheatwaves.constants as const
import ehfheatwaves.getoptions as getoptions
import ehfheatwaves.ncio as ncio
import ehfheatwaves.qtiler as qtiler
from ehfheatwaves.getoptions import options

warnings.simplefilter('ignore', category=RuntimeWarning)


class GridDescription(object):
    """Description of grid."""

    def __init__(self, lats=np.array([])):
        """## Arguments:

        - lats : Latitudes.
        """
        self.lats = lats


def window_percentile(temp:np.ndarray, daysinyear:int=365, wsize:int=15)->np.ndarray:
    """Calculate a day-of-year moving window percentile.

    ## Arguments:

    - temp : Temperature data.
    - daysinyear : Number of days in a year.
    - wsize : Number of days in moving window.
    """
    # Initialise array.
    pctl = np.ones(((daysinyear,)+temp.shape[1:]))*const.FILL_VAL

    # Construct the window.
    window = np.zeros(daysinyear, dtype=bool)
    window[-np.floor(wsize/2.).astype(int):] = 1
    window[:np.ceil(wsize/2.).astype(int)] = 1
    window = np.tile(window, options.bpend + 1 - options.bpstart)

    # Select the interpolation method.
    if options.qtilemethod=='python':
        percentile = np.percentile
        parameter = 0
    elif options.qtilemethod=='zhang':
        percentile = qtiler.quantile_zhang_fast
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
        window = np.roll(window, 1)

    # Remaining nans are missing data.
    pctl[np.isnan(pctl)] = const.MISSING_VAL

    return pctl


def identify_hw(ehfs:np.ndarray)->tuple:
    """Locate heatwaves from EHF and returns an event indicator and a duration indicator.

    ## Arguments:

    - ehfs : EHF values.

    ## Returns:
    - events : array of bools for heatwave events.
    - endss : array of integers for heatwave duration.
    """
    # Agregate consecutive days with EHF>0
    # First day contains duration
    if np.isnan(ehfs).any(): # then ehfs is a view to tmax or tmin
        # This is handled differently to EHFs because tmax or tmin could theoretically hold -ve
        # values. Therefore we identify values of the exceedence index that are not nan values as
        # heatwaves.
        events = np.logical_not(np.isnan(ehfs)).astype(int)
    else: # ehfs is a view to actual EHF values
        events = (ehfs>0.0).astype(int)
        events[events.mask==True] = 0
    for i in range(events.shape[0] - 2, -1, -1):
         events[i,events[i,...]>0] = events[i+1,events[i,...]>0]+1

    # Identify when heatwaves start with duration
    # Given that first day contains duration
    diff = np.zeros(events.shape)
    # Insert the first diff value as np.diff doesn't catch it because
    # there is no pevious value to compare to.
    diff[0,...] = events[0,...]
    diff[1:,...] = np.diff(events, axis=0)
    endss = np.ma.zeros(ehfs.shape, dtype=int)
    endss[diff>2] = events[diff>2]

    # Remove events less than 3 days
    events[diff==2] = 0
    events[np.roll(diff==2, 1, axis=0)] = 0
    events[diff==1] = 0
    del diff
    events[events>0] = 1
    events = events.astype(bool)
    endss[endss<3] = 0
    return events, endss


def identify_semi_hw(ehfs:np.ndarray)->tuple:
    """identify_hw locates heatwaves from EHF and returns an event indicator and a duration
    indicator. This function does not exclude events less than three days in duration.

    ## Arguments:

    - ehfs : EHF values.

    ## Returns:

    - events : array of bools for heatwave events.
    - endss -- array of integers for heatwave duration.
    """
    # Agregate consecutive days with EHF>0
    # First day contains duration
    if np.isnan(ehfs).any(): # then ehfs is a view to tmax or tmin
        # This is handled differently to EHFs because tmax or tmin could theoretically hold -ve
        # values. Therefore we identify values of the exceedence index that are not nan values as
        # heatwaves.
        events = np.logical_not(np.isnan(ehfs)).astype(int)
    else: # ehfs is a view to actual EHF values
        events = (ehfs>0.0).astype(int)
        events[events.mask==True] = 0
    for i in range(events.shape[0] - 2, -1, -1):
         events[i,events[i,...]>0] = events[i+1,events[i,...]>0]+1

    # Identify when heatwaves start with duration
    # Given that first day contains duration
    diff = np.zeros(events.shape)
    # Insert the first diff value as np.diff doesn't catch it because
    # there is no pevious value to compare to.
    diff[0,...] = events[0,...]
    diff[1:,...] = np.diff(events, axis=0)
    endss = np.ma.zeros(ehfs.shape, dtype=int)
    endss[diff>0] = events[diff>0]
    del diff
    events[events>0] = 1
    events = events.astype(bool)
    return events, endss


def hw_aspects(EHF:np.ndarray, season:str, hemisphere:str)->tuple:
    """Call `identify_hw` and/or `identify_semi_hw` and calculate seasonal aspects.

    ## Arguments:

    - EHF : EHF values.
    - season : The season in which to calculate aspects for.
    - hemisphere : The hemisphere to calculate heatwaves for.

    ## Returns:

    - HWA : Amplitude
    - HWM : Magnitude
    - HWN : Number
    - HWF : Frequency
    - HWD : Maximum duration
    - HWT : Timing
    """
    global timedata
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
    HWA = np.ones(((timedata.nyears,)+(EHF.shape[1],)))*const.FILL_VAL
    HWM = HWA.copy()
    HWN = HWA.copy()
    HWF = HWA.copy()
    HWD = HWA.copy()
    HWT = HWA.copy()
    # Loop over years
    for iyear, year in enumerate(range(timedata.first_year, timedata.daylast.year)):
        if options.oldmethod:
            if (year==timedata.daylast.year): continue # Incomplete yr
            # Select this years season
            allowance = 14 # For including heawave days after the end of the season
            ifrom = startday + timedata.daysinyear*iyear - 1 # -1 to include Oct 31st
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
        HWT[iyear,...] = np.argmax(event_i, axis=0)
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
        # Locate invalid values or misisng values
        missing = EHF_i.mask.all(axis=0)
        if missing.any():
            HWT[iyear,missing] = const.MISSING_VAL
            HWN[iyear,missing] = const.MISSING_VAL
            HWF[iyear,missing] = const.MISSING_VAL
            HWD[iyear,missing] = const.MISSING_VAL
            HWA[iyear,missing] = const.MISSING_VAL
            HWM[iyear,missing] = const.MISSING_VAL
        invalid = HWN[iyear,...]==0
        HWT[iyear,invalid] = const.INVALID_VAL
        HWD[iyear,invalid] = const.INVALID_VAL
        HWA[iyear,invalid] = const.INVALID_VAL
        HWM[iyear,invalid] = const.INVALID_VAL
    return HWA, HWM, HWN, HWF, HWD, HWT


# Calculate metrics year by year
def split_hemispheres(EHF:np.ndarray, north:bool, south:bool)->tuple:
    """Split the input data by hemispheres, and glue them back together after heatwave
    calculations.

    The EHF spatial axes are reshaped into a single dimension.
    The output arrays are 2D. When saving, data should be reshaped or indexed with a land-sea mask.

    ## Arguments:

    - EHF : EHF values.
    - north : Calculate heatwaves for the northern hemisphere.
    - south : Calculate heatwaves for the southern hemisphere.

    ## Returns:

    - HWA : Amplitude
    - HWM : Magnitude
    - HWN : Number
    - HWF : Frequency
    - HWD : Maximum duration
    - HWT : Timing
    """
    lats = grid.lats
    if south:
        if options.maskfile:
            EHF_s = EHF[:,:(mask[lats<=0]>0).sum()]
        else:
            EHF_s = EHF[:,lats<=0,...]
        # Reshape to 2D
        space = EHF_s.shape[1:]
        if len(space)>1:
            EHF_s = EHF_s.reshape(EHF_s.shape[0], space[0]*space[1])
        # Southern hemisphere aspects
        HWA_s, HWM_s, HWN_s, HWF_s, HWD_s, HWT_s = hw_aspects(EHF_s, options.season, 'south')
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
        HWA_n, HWM_n, HWN_n, HWF_n, HWD_n, HWT_n = hw_aspects(EHF_n, options.season, 'north')
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


def main():
    """Main function that is called by the entry point script.

    See the following documentation for more information on entry points.
    https://setuptools.pypa.io/en/latest/userguide/entry_point.html?highlight=entry

    This function parses the command line arguments and options then calculates the heatwaves and
    saves them to netcdf files.
    """
    global grid, timedata, mask, options
    # Get the options and variables
    options = getoptions.parse_arguments(sys.argv[1:])

    # Load time data
    if options.verbose: print("Loading data")
    if options.tmaxfile:
        filename = options.tmaxfile
    elif options.tminfile:
        filename = options.tminfile
    timedata = ncio.TimeData(filename, options.timevname)
    grid = GridDescription()

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
        tpct = window_percentile(tave_base, daysinyear=timedata.daysinyear)
    if options.keeptmax:
        txpct = window_percentile(tmax, daysinyear=timedata.daysinyear)
    if options.keeptmin:
        tnpct = window_percentile(tmin, daysinyear=timedata.daysinyear)
    if not options.noehf: del tave_base

    # Load all data
    if options.verbose: print("Loading data")
    if options.keeptave or options.keeptmax:
        tmax, grid.lats = ncio.get_all_data(options.tmaxfile, options.tmaxvname, options)
        original_shape = tmax.shape
    if options.keeptave or options.keeptmin:
        tmin, grid.lats = ncio.get_all_data(options.tminfile, options.tminvname, options)
        original_shape = tmin.shape
    if options.keeptave:
        tave = (tmax + tmin)/2.

    # Remove leap days from data
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
        timedata.calendar = '365_day'

    # Remove incomplete starting year
    timedata.first_year = timedata.dayone.year
    if (timedata.dayone.month!=1)|(timedata.dayone.day!=1):
        timedata.first_year = timedata.dayone.year+1
        start = np.argmax(timedata.dates.year==timedata.first_year)
        if options.keeptave:
            tave = tave[start:,...]
            original_shape = (tave.shape[0], original_shape[1], original_shape[2])
        if options.keeptmax:
            tmax = tmax[start:,...]
            original_shape = (tmax.shape[0], original_shape[1], original_shape[2])
        if options.keeptmin:
            tmin = tmin[start:,...]
            original_shape = (tmin.shape[0], original_shape[1], original_shape[2])

    # Apply mask
    if options.maskfile:
        if options.keeptave:
            tave = tave[:,mask]
        if options.keeptmax:
            tmax = tmax[:,mask]
        if options.keeptmin:
            tmin = tmin[:,mask]

    if options.verbose: print("Caclulating definition")
    # Calculate EHF
    if not options.noehf:
        EHF = np.ma.ones(tave.shape)*np.nan
        for i in range(32, tave.shape[0]):
            EHIaccl = tave[i-2:i+1,...].sum(axis=0)/3.0 - tave[i-32:i-2,...].sum(axis=0)/30.0
            EHIsig = tave[i-2:i+1,...].sum(axis=0)/3.0 - tpct[i-timedata.daysinyear*int((i+1)/timedata.daysinyear),...]
            EHF[i,...] = np.ma.maximum(EHIaccl, 1.0)*EHIsig
        EHF[EHF<0] = 0
    if options.ehi:
        EHIaccl = np.ma.ones(tave.shape)*np.nan
        EHIsig = np.ma.ones(tave.shape)*np.nan
        for i in range(32, tave.shape[0]):
            EHIaccl[i,...] = tave[i-2:i+1,...].sum(axis=0)/3.0 - tave[i-32:i-2,...].sum(axis=0)/30.0
            EHIsig[i,...] = tave[i-2:i+1,...].sum(axis=0)/3.0 - tpct[i-timedata.daysinyear*int((i+1)/timedata.daysinyear),...]

    # Tx90pc exceedences
    if options.keeptmin or options.keeptmax:
        if options.keeptmax:
            txexceed = np.ma.ones(tmax.shape)*np.nan
            for i in range(0, tmax.shape[0]):
                idoy = i-timedata.daysinyear*int((i+1)/timedata.daysinyear)
                txexceed[i,...] = tmax[i,...]>txpct[idoy,...]
            txexceed[txexceed>0] = tmax[txexceed>0]
            txexceed[txexceed==0] = np.nan # need nans to be identified correctly
        if options.keeptmin:
            tnexceed = np.ma.ones(tmin.shape)*np.nan
            for i in range(0, tmin.shape[0]):
                idoy = i-timedata.daysinyear*int((i+1)/timedata.daysinyear)
                tnexceed[i,...] = tmin[i,...]>tnpct[idoy,...]
            tnexceed[tnexceed>0] = tmin[tnexceed>0]
            tnexceed[tnexceed==0] = np.nan # need nans to be identified correctly

    # Calculate daily output
    if options.dailyout or options.dailyout: event, ends = identify_hw(EHF)
    if options.tx90pcd: event_tx, ends_tx = identify_hw(txexceed)
    if options.tn90pcd: event_tn, ends_tn = identify_hw(tnexceed)

    timedata.nyears = len(range(timedata.first_year, timedata.daylast.year + 1))

    # Calculate yearly output
    if options.yearlyout:
        if options.verbose: print("Calculating yearly aspects")
        # Split by latitude
        north = (grid.lats>0).any()
        south = (grid.lats<=0).any()
        if not options.noehf:
            HWA_EHF, HWM_EHF, HWN_EHF, HWF_EHF, HWD_EHF, HWT_EHF = split_hemispheres(
                    EHF,
                    north,
                    south,
                    )
        if options.tx90pc:
            HWA_tx, HWM_tx, HWN_tx, HWF_tx, HWD_tx, HWT_tx = split_hemispheres(
                    txexceed,
                    north,
                    south,
                    )
        if options.tn90pc:
            HWA_tn, HWM_tn, HWN_tn, HWF_tn, HWD_tn, HWT_tn = split_hemispheres(
                    tnexceed,
                    north,
                    south,
                    )

    if options.verbose: print("Saving")
    # Save yearly data to netcdf
    if options.yearlyout:
        if not options.noehf:
            ncio.save_yearly(
                    HWA_EHF,
                    HWM_EHF,
                    HWN_EHF,
                    HWF_EHF,
                    HWD_EHF,
                    HWT_EHF,
                    tpct,
                    "EHF",
                    timedata,
                    options,
                    mask,
                    )
        if options.tx90pc:
            ncio.save_yearly(
                    HWA_tx,
                    HWM_tx,
                    HWN_tx,
                    HWF_tx,
                    HWD_tx,
                    HWT_tx,
                    txpct,
                    "tx90pct",
                    timedata,
                    options,
                    mask,
                    )
        if options.tn90pc:
            ncio.save_yearly(
                    HWA_tn,
                    HWM_tn,
                    HWN_tn,
                    HWF_tn,
                    HWD_tn,
                    HWT_tn,
                    tnpct,
                    "tn90pct",
                    timedata,
                    options,
                    mask,
                    )

    # Save daily data to netcdf
    if options.dailyout:
        if options.keeptave:
            ncio.save_daily(
                    EHF,
                    event,
                    ends,
                    options,
                    timedata,
                    original_shape,
                    mask,
                    defn='EHF',
                    )
        if options.tx90pcd:
            ncio.save_daily(
                    txexceed,
                    event_tx,
                    ends_tx,
                    options,
                    timedata,
                    original_shape,
                    mask,
                    defn='tx90pct',
                    )
        if options.tn90pcd:
            ncio.save_daily(
                    tnexceed,
                    event_tn,
                    ends_tn,
                    options,
                    timedata,
                    original_shape,
                    mask,
                    defn='tn90pct',
                    )

    # save EHIs
    if options.ehi:
        ncio.save_ehi(EHIsig, EHIaccl, options, timedata, original_shape, mask)

