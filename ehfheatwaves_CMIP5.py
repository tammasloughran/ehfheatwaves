import warnings
warnings.filterwarnings("ignore")
import sys
try: 
    import pandas as pd
except ImportError:
    print "Please install pandas"
    sys.exit(2)
import numpy as np
import datetime as dt
import math
import qtiler
from netCDF4 import MFDataset, Dataset
import netcdftime
from optparse import OptionParser
import pdb

# Parse command line arguments
usage = 'usage: %prog -x <FILE> -n <FILE> [-a <FILE>] -m <FILE> [options]'
parser = OptionParser(usage=usage)
parser.add_option('-x', '--tmax', dest='tmaxfile', 
        help='file containing tmax', metavar='FILE')
parser.add_option('-n', '--tmin', dest='tminfile',
        help='file containing tmin', metavar='FILE')
parser.add_option('-m', '--mask', dest='maskfile',
        help='file containing land-sea mask', metavar='FILE')
parser.add_option('-s', '--season', dest='season', default='summer',
        help='austal season for annual metrics. Defaults to austral summer',
        metavar='STR')
parser.add_option('-p', dest='pcntl', type='float', default=90,
        help='the percentile to use for thresholds. Defaults to 90',
        metavar='INT')
parser.add_option('--base', dest='bp', default='1961-1990',
        help='base period to calculate thresholds. Default 1961-1990',
        metavar='YYYY-YYYY')
parser.add_option('-q', '--qmethod', dest='qtilemethod', default='climpact',
        help='quantile interpolation method. Default is climpact', 
        metavar='STR')
parser.add_option('-d', '--daily', action="store_true", dest='daily', 
        help='output daily EHF values and heatwave indicators')
parser.add_option('--dailyonly', action="store_true", dest='dailyonly',
        help='output only daily values and suppress yearly output')
(options, args) = parser.parse_args()
if not options.tmaxfile or not options.tminfile:
    print 'Please specify tmax and tmin files if not using tave.'
    sys.exit(2)
if not options.maskfile:
    print 'Please specify a land-sea mask file.'
    sys.exit(2)
if len(options.bp)!=9:
    print 'Incorect base period format.'
    sys.exit(2)
else:
    bpstart = int(options.bp[:4])
    bpend = int(options.bp[5:9])
# Percentile
pcntl = options.pcntl
# climpact/python/matlab
qtilemethod = options.qtilemethod
# season (winter/summer)
season = options.season
if (season!='summer')&(season!='winter'):
    print 'Use either summer or winter. (Austral)'
    sys.exit(2)
# save daily EHF output
dailyout = options.daily
if options.dailyonly: 
    dailyout = True
    yearlyout = False

# Load data
tmaxnc = MFDataset(options.tmaxfile, 'r')
nctime = tmaxnc.variables['time']
calendar = nctime.calendar
if (calendar=='gregorian')|(calendar=='proleptic_gregorian')|(calendar=='365_day'):
    daysinyear = 365
    if season=='winter':
        seasonlen = 153
        startday = 121
        endday = 274
    else:
        seasonlen = 151
        startday = 304
        endday = 455
    dayone = netcdftime.num2date(nctime[0], nctime.units, 
            calendar=calendar)
    daylast = netcdftime.num2date(nctime[-1], nctime.units, 
            calendar=calendar)
    pdb.set_trace()
    dates = pd.date_range(dayone, daylast)
elif calendar=='360_day':
    daysinyear = 360
    seasonlen = 150
    if season=='winter':
        startday = 121
        endday = 271
    else:
        startday = 301
        endday = 451
    dayone = netcdftime.num2date(nctime[0], nctime.units,
            calendar=calendar)
    daylast = netcdftime.num2date(nctime[-1], nctime.units,
            calendar=calendar)
    class calendar360():
        def __init__(self,sdate,edate):
            self.year = np.repeat(range(sdate.year,edate.year+1), 360, 0)
            self.month = np.tile(np.repeat(range(1,12+1), 30, 0), 104)
            self.day = np.tile(np.tile(range(1,30+1), 12), 104)
            sdoyi = (sdate.month-1)*30+sdate.day
            self.year = self.year[sdoyi:]
            self.month = self.month[sdoyi:]
            self.day = self.day[sdoyi:]
            edoyi = (12-edate.month)*30+(31-edate.day)
            self.year = self.year[:-edoyi]
            self.month = self.month[:-edoyi]
            self.day = self.day[:-edoyi]
    dates = calendar360(dayone, daylast)
else:
    print 'Unrecognized calendar'
    sys.exit(2)

# Load land-sea mask
masknc = Dataset(options.maskfile, 'r')
mask = masknc.variables['sftlf'][:]
mask = mask.astype(np.bool)
masknc.close()

# Load base period data
tmax = tmaxnc.variables['tasmax'][(bpstart<=dates.year)&(dates.year<=bpend)]
original_shape = tmax.shape
tmax = tmax[:,mask]
if tmaxnc.variables['tasmax'].units=='K': tmax -= 273.15
tminnc = MFDataset(options.tminfile, 'r')
tmin = tminnc.variables['tasmin'][(bpstart<=dates.year)&(dates.year<=bpend)]
tmin = tmin[:,mask]
if tminnc.variables['tasmin'].units=='K': tmin -= 273.15
tave_base = (tmax + tmin)/2.
del tmin
del tmax

# Remove leap days in gregorian calendars
if (calendar=='gregorian')|(calendar=='proleptic_gregorian'):
    dates_base = dates[(bpstart<=dates.year)&(dates.year<=bpend)]
    tave_base = tave_base[(dates_base.month!=2)|(dates_base.day!=29),...]
    del dates_base

# Caclulate 90th percentile
tpct = np.ones((daysinyear,tave_base.shape[1]))*np.nan
window = np.zeros(daysinyear,dtype=np.bool)
wsize = 15.
window[-np.floor(wsize/2.):] = 1
window[:np.ceil(wsize/2.)] = 1
window = np.tile(window,30)
if qtilemethod=='python':
    percentile = np.percentile
    parameter = 0
elif qtilemethod=='zhang':
    percentile = qtiler.quantile_zhang
    parameter = False
elif qtilemethod=='matlab':
    percentile = qtiler.quantile_R
    parameter = 5
elif qtilemethod=='climpact':
    percentile = qtiler.quantile_climpact
    parameter = False
for day in xrange(daysinyear):
    tpct[day,...] = percentile(tave_base[window,...], pcntl, parameter)
    window = np.roll(window,1)
del tave_base
del window

# Load data
tmax = tmaxnc.variables['tasmax'][:]
tmin = tminnc.variables['tasmin'][:]
tmax = tmax[:,mask]
if tmaxnc.variables['tasmax'].units=='K': tmax -= 273.15
tmin = tmin[:,mask]
if tminnc.variables['tasmin'].units=='K': tmin -= 273.15
tave = (tmax + tmin)/2.
del tmax
del tmin

# Remove incomplete starting year
if (dayone.month!=1)|(dayone.day!=1):
    i = np.argmax(dates.year==dayone.year+1)
    tave = tave[i:,...]

# Calculate EHF
EHF = np.ones(tave.shape)*np.nan
for i in xrange(32,tave.shape[0]):
    EHIaccl = tave[i-2:i+1,...].sum(axis=0)/3. - \
            tave[i-32:i-2,...].sum(axis=0)/30.
    EHIsig = tave[i-2:i+1,...].sum(axis=0)/3. - \
            tpct[i-365*(i/364),...]
    EHF[i,...] = np.maximum(EHIaccl,1.)*EHIsig
EHF[EHF<0] = 0

# Agregate consecutive days with EHF>0
# First day contains duration
event = (EHF>0.).astype(np.int)
for i in xrange(event.shape[0]-2,0,-1):
    event[i,event[i,...]>0] = event[i+1,event[i,...]>0]+1

# Identify when heatwaves start with duration
# Given that first day contains duration
diff = np.zeros(event.shape)
diff[1:,...] = np.diff(event, axis=0)
ends = np.zeros(tave.shape,dtype=np.int)
ends[diff>2] = event[diff>2]
del diff

# Remove events less than 3 days
event[ends==2] = 0
event[np.roll(ends==2, 1, axis=0)] = 0
event[ends==1] = 0
event[event>0] = 1
event = event.astype(np.bool)
ends[ends<3] = 0

# Calculate metrics year by year
nyears = len(range(dayone.year,daylast.year))
HWN = np.ones((nyears,tave.shape[1]))*np.nan
HWF = HWN.copy()
HWD = HWN.copy()
HWA = HWN.copy()
HWM = HWN.copy()
HWT = HWN.copy()
for iyear, year in enumerate(xrange(dayone.year,daylast.year)):
    if (year==daylast.year)&(season=='summer'): continue # Incomplete last yr
    HWMtmp = np.ones((daysinyear,tave.shape[1]))*np.nan
    # Select this years season
    ifrom = startday+daysinyear*iyear
    ito = endday+daysinyear*iyear
    EHF_i = EHF[ifrom:ito,...]
    event_i = event[ifrom:ito,...]
    duration_i = ends[ifrom:ito,...].astype(np.float)
    # Calculate metrics
    HWN[iyear,...] = (duration_i>0).sum(axis=0)
    HWF[iyear,...] = 100*duration_i.sum(axis=0)/float(seasonlen)
    HWD[iyear,...] = duration_i.max(axis=0)
    HWA[iyear,...] = EHF_i.max(axis=0)
    HWMtmp[event_i==True] = EHF_i[event_i==True]
    HWM[iyear,...] = np.nanmean(HWMtmp, axis=0)
    HWT[iyear,...] = np.argmax(event_i,axis=0)

# Save to netCDF
if yearlyout:
    experiment = tmaxnc.__getattribute__('experiment')
    model = tmaxnc.__getattribute__('model_id')
    realization = tmaxnc.__getattribute__('realization')
    yearlyout = Dataset('EHF_heatwaves_%s_%s_%s_yearly_%s.nc'%(model, 
            experiment, realization, season), mode='w')
    yearlyout.createDimension('time', len(range(dayone.year,
            daylast.year)))
    yearlyout.createDimension('lon', tmaxnc.dimensions['lon'].__len__())
    yearlyout.createDimension('lat', tmaxnc.dimensions['lat'].__len__())
    yearlyout.createDimension('day', daysinyear)
    setattr(yearlyout, "author", "Tammas Loughran")
    setattr(yearlyout, "contact", "t.loughran@student.unsw.edu.au")
    setattr(yearlyout, "date", dt.datetime.today().strftime('%Y-%m-%d'))
    setattr(yearlyout, "script", "ehfheatwaves_CMIP5.py")
    setattr(yearlyout, "model_id", model)
    setattr(yearlyout, "experiment", experiment)
    setattr(yearlyout, "base_period", "%s-%s"%(str(bpstart),str(bpend)))
    setattr(yearlyout, "percentile", "%sth"%(str(pcntl)))
    setattr(yearlyout, "realization", "%s"%(realization))
    otime = yearlyout.createVariable('time', 'f8', 'time', 
            fill_value=-999.99)
    setattr(otime, 'units', 'year')
    olat = yearlyout.createVariable('lat', 'f8', 'lat')
    setattr(olat, 'Longname', 'Latitude')
    setattr(olat, 'units', 'degrees_north')
    olon = yearlyout.createVariable('lon', 'f8', 'lon')
    setattr(olon, 'Longname', 'Longitude')
    setattr(olon, 'units', 'degrees_east')
    otpct = yearlyout.createVariable('t%spct'%(pcntl), 'f8', 
	    ('day','lat','lon'), fill_value=-999.99)
    setattr(otpct, 'Longname', '90th percentile')
    setattr(otpct, 'units', 'degC')
    setattr(otpct, 'description', 
            '90th percentile of %s-%s'%(str(bpstart),str(bpend)))
    HWAout = yearlyout.createVariable('HWA_EHF', 'f8', ('time','lat','lon'), 
            fill_value=-999.99)
    setattr(HWAout, 'Longname', 'Peak of the hottest heatwave per year')
    setattr(HWAout, 'units', 'degC2')
    setattr(HWAout, 'description', 
            'Peak of the hottest heatwave per year')
    HWMout = yearlyout.createVariable('HWM_EHF', 'f8', ('time','lat','lon'),
            fill_value=-999.99)
    setattr(HWMout, 'Longname', 'Average magnitude of the yearly heatwave')
    setattr(HWMout, 'units', 'degC2')
    setattr(HWMout, 'description', 'Average magnitude of the yearly heatwave')
    HWNout = yearlyout.createVariable('HWN_EHF', 'f8', ('time', 'lat', 'lon'), 
            fill_value=-999.99)
    setattr(HWNout, 'Longname', 'Number of heatwaves')
    setattr(HWNout, 'units','')
    setattr(HWNout, 'description', 'Number of heatwaves per year')
    HWFout = yearlyout.createVariable('HWF_EHF', 'f8', ('time','lat','lon'), 
            fill_value=-999.99)
    setattr(HWFout,'Longname','Number of heatwave days')
    setattr(HWFout, 'units', '%')
    setattr(HWFout, 'description', 'Proportion of heatwave days per season')
    HWDout = yearlyout.createVariable('HWD_EHF', 'f8', ('time','lat','lon'), 
            fill_value=-999.99)
    setattr(HWDout, 'Longname', 'Duration of yearly longest heatwave')
    setattr(HWDout, 'units', 'days')
    setattr(HWDout, 'description', 'Duration of the longest heatwave per year')
    HWTout = yearlyout.createVariable('HWT_EHF', 'f8', ('time','lat','lon'), 
            fill_value=-999.99)
    setattr(HWTout, 'Longname', 'First heat wave day of the season')
    otime[:] = range(dayone.year, daylast.year)
    olat[:] = tmaxnc.variables['lat'][:]
    olon[:] = tmaxnc.variables['lon'][:]
    dummy_array = np.ones((daysinyear,)+original_shape[1:])*np.nan
    dummy_array[:,mask] = tpct
    otpct[:] = dummy_array.copy()
    dummy_array = np.ones((nyears,)+original_shape[1:])*np.nan
    dummy_array[:,mask] = HWA
    HWAout[:] = dummy_array.copy()
    dummy_array[:,mask] = HWM
    HWMout[:] = dummy_array.copy()
    dummy_array[:,mask] = HWN
    HWNout[:] = dummy_array.copy()
    dummy_array[:,mask] = HWF
    HWFout[:] = dummy_array.copy()
    dummy_array[:,mask] = HWD
    HWDout[:] = dummy_array.copy() 
    dummy_array[:,mask] = HWT
    HWTout[:] = dummy_array.copy()
    yearlyout.close()

if dailyout:
    dailyout = Dataset('EHF_heatwaves_%s_%s_%s_daily.nc'%(model, experiment, realization), mode='w')
    dailyout.createDimension('time', tmaxnc.dimensions['time'].__len__())
    dailyout.createDimension('lon', tmaxnc.dimensions['lon'].__len__())
    dailyout.createDimension('lat', tmaxnc.dimensions['lat'].__len__())
    setattr(dailyout, "author", "Tammas Loughran")
    setattr(dailyout, "contact", "t.loughran@student.unsw.edu.au")
    setattr(dailyout, "date", dt.datetime.today().strftime('%Y-%m-%d'))
    setattr(dailyout, "script", "ehfheatwaves_CMIP5.py")
    setattr(dailyout, "base_period", "%s-%s"%(str(bpstart),str(bpend)))
    setattr(dailyout, "percentile", "%sth"%(str(pcntl)))
    setattr(dailyout, "model_id", model)
    setattr(dailyout, "experiment", experiment)
    setattr(dailyout, "realization", realization)
    otime = dailyout.createVariable('time', 'f8', 'time',
                    fill_value=-999.99)
    setattr(otime, 'units', tmaxnc.variables['time'].units)
    setattr(otime, 'calendar', tmaxnc.variables['time'].calendar)
    olat = dailyout.createVariable('lat', 'f8', 'lat')
    setattr(olat, 'Longname', 'Latitude')
    setattr(olat, 'units', 'degrees_north') 
    olon = dailyout.createVariable('lon', 'f8', 'lon')
    setattr(olon, 'Longname', 'Longitude')
    setattr(olon, 'units', 'degrees_east')
    oehf = dailyout.createVariable('ehf', 'f8', ('time','lat','lon'),
                fill_value=-999.99)
    setattr(oehf, 'Longname', 'Excess Heat Factor')
    setattr(oehf, 'units', 'degC2')
    oevent = dailyout.createVariable('event', 'f8', ('time','lat','lon'),
                fill_value=-999.99)
    setattr(oehf, 'Longname', 'Event indicator')
    oends = dailyout.createVariable('ends', 'f8', ('time','lat','lon'),
                        fill_value=-999.99)
    setattr(oends, 'Longname', 'Duration at start of heatwave')
    setattr(oends, 'units', 'days')
    otime[:] = tmaxnc.variables['time'][:]
    olat[:] = tmaxnc.variables['lat'][:]
    olon[:] = tmaxnc.variables['lon'][:]
    dummy_array = np.ones((EHF.shape[0],)+original_shape[1:])*np.nan
    dummy_array[:,mask] = EHF
    oehf[:] = dummy_array.copy()
    dummy_array[:,mask] = event
    oevent[:] = dummy_array.copy()
    dummy_array[:,mask] = ends
    oends[:] = dummy_array.copy()
    dailyout.close()
