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
parser.add_option('-a', '--tave', dest='tavefile',
        help='file containing tave', metavar='FILE')
parser.add_option('-m', '--mask', dest='maskfile',
        help='file containing land-sea mask', metavar='FILE')
parser.add_option('-s', '--season', dest='season', default='summer',
        help='season for annual metrics. Defaults to austral summer',
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
parser.add_option('-y', '--yearly', action="store_true", dest='yearly',
        help='output yearly heatwaaspects')
(options, args) = parser.parse_args()
if options.tavefile:
    usetave = True
elif not options.tmaxfile or not options.tminfile:
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
# save daily EHF output
dailyout = options.daily
yearlyout = options.yearly

# Load data
tmaxnc = MFDataset(options.tmaxfile, 'r')
nctime = tmaxnc.variables['time']
calendar = nctime.calendar
if calendar=='gregorian' or 'proleptic_gregorian' or '365_day':
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
    class calendar360:
        def __init__(self,sdate,edate):
            self.year = np.repeat(range(sdate.year,edate.year+1), 360, 0)
            self.month = np.tile(np.repeat(range(1,12+1), 30, 0), 104)
            self.day = np.tile(np.tile(range(1,30+1), 12), 104)
            sdoyi = (sdate.month-1)*30+sdate.day
            self.year = self.year[sdoy:]
            self.month = self.month[sdoy:]
            self.day = self.day[sdoy:]
            edoyi = (12-edate.month)*30+(31-edate.day)
            self.year = self.year[:-edoyi]
            self.month = self.month[:-edoyi]
            self.day = self.day[:-edoyi]
    dates = calendar360(dayone, daylast)

# Load land-sea mask
masknc = Dataset(options.maskfile, 'r')
mask = masknc.variables['sftlf'][:]
mask = mask.astype(np.bool)
masknc.close()

# Load base period data
tmax = tmaxnc.variables['tasmax'][(bpstart<=dates.year)&(dates.year<=bpend)]
original_shape = tmax.shape
tmax = tmax[:,mask]
tminnc = MFDataset(options.tminfile)
tmin = tminnc.variables['tasmin'][(bpstart<=dates.year)&(dates.year<=bpend)]
tmin = tmin[:,mask]
tave_base = (tmax + tmin)/2.
del tmin
del tmax

# Remove leap days in gregorian calendars
if calendar=='gregorian' or 'proleptic_gregorian':
    dates_base = dates[(bpstart<=dates.year)&(dates.year<=bpend)]
    tave_base = tave_base[(dates_base.month!=2)|(dates_base.day!=29),...]
    del dates_base

# Select Base period
#base = ((bpstart<=dates.year)&(dates.year<=bpend))
#tave_base = tave[base,...]

print 'calculating %ile'
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

#pdb.set_trace()

# Create a netcdf file for EHF
dailyout = Dataset('EHF_heatwaves_ACCESS_daily.nc', mode='w')
dailyout.createDimension('time', len(dates))
dailyout.createDimension('lon', tmaxnc.dimensions['lon'].__len__())
dailyout.createDimension('lat', tmaxnc.dimensions['lat'].__len__())
setattr(dailyout, "author", "Tammas Loughran")
setattr(dailyout, "contact", "t.loughran@student.unsw.edu.au")
setattr(dailyout, "date", dt.datetime.today().strftime('%Y-%m-%d'))
setattr(dailyout, "script", "ehfheatwaves.py")
setattr(dailyout, "base_period", "%s-%s"%(str(bpstart),str(bpend)))
setattr(dailyout, "percentile", "%sth"%(str(pcntl)))
setattr(dailyout, "dataset", "ACCESS")
otime = dailyout.createVariable('time', 'f8', 'time',
                fill_value=-999.99)
setattr(otime, 'units', 'day as %Y%m%d.%f')
setattr(otime, 'calendar', '365_day')
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

daysinchunks = 50.*daysinyear
nchunks = int(math.floor(len(dates)/daysinchunks))
finalchunki = int(nchunks*daysinchunks)

# Load 50 year chunks
for chunk in xrange(nchunks):
    print 'loading chunk', chunk
    i1 = int(daysinchunks*chunk)
    i2 = int(daysinchunks*(chunk+1))
    if chunk==0:
        tmax = tmaxnc.variables['tasmax'][i1:i2]
        tmin = tminnc.variables['tasmin'][i1:i2]
    elif chunk==nchunks:
        tmax = tmaxnc.variables['tasmax'][finalchunki-31:]
        tmin = tminnc.variables['tasmin'][finalchunki-31:]
    else:
        tmax = tmaxnc.variables['tasmax'][i1-31:i2]
        tmin = tminnc.variables['tasmin'][i1-31:i2]
    tmax = tmax[:,mask]
    tmin = tmin[:,mask]
    tave = (tmax + tmin)/2.
    del tmax
    del tmin
    print 'calculating ehf'
    # Calculate EHF
    EHF = np.ones(tave.shape)*np.nan
    for i in xrange(32,tave.shape[0]):
        #if i<31: continue
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
    # Last day contains duration.
    #for i in range(1,event.shape[0]):
    #    event[i,event[i,...]>0] = event[i-1,event[i,...]>0]+1
    print 'identifying duration'
    # Identify when heatwaves start with duration
    # Given that first day contains duration
    diff = np.zeros(event.shape)
    diff[1:,...] = np.diff(event, axis=0)
    ends = np.zeros(tave.shape,dtype=np.int)
    ends[diff>2] = event[diff>2]
    del diff
    print 'removing junk events'
    # Remove events less than 3 days
    event[ends==2] = 0
    event[np.roll(ends==2, 1, axis=0)] = 0
    event[ends==1] = 0
    event[event>0] = 1
    event = event.astype(np.bool)
    ends[ends<3] = 0
    print 'saving'
    # Save EHF to file
    #dummy_array = np.ones((EHF.shape[0],)+original_shape[1:])*np.nan
    dummy_array = np.ones(original_shape[1:])*np.nan
    for i in xrange(32, EHF.shape[0]):
        dummy_array[mask] = EHF[i,...]
        oehf[i+i1] = dummy_array.copy()
        del EHF
        dummy_array[mask] = event
        oevent[i+i1] = dummy_array.copy()
        del event
        dummy_array[mask] = ends
        oends[i+i1] = dummy_array.copy()
        del ends

dailyout.close()
sys.exit()



# Calculate metrics year by year
nyears = len(range(syear,eyear))
HWN = np.ones((nyears,tave.shape[1]))*np.nan
HWF = HWN.copy()
HWD = HWN.copy()
HWA = HWN.copy()
HWM = HWN.copy()
HWT = HWN.copy()
for iyear, year in enumerate(xrange(syear,eyear)):
    if (year==eyear)&(season=='summer'): continue # Last year is incomplete
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
    yearlyout = Dataset('EHF_heatwaves_1911-2014_yearly_%s.nc'%(season), mode='w')
    yearlyout.createDimension('time', len(range(syear,eyear)))
    yearlyout.createDimension('lon', len(lons))
    yearlyout.createDimension('lat', len(lats))
    yearlyout.createDimension('day', daysinyear)
    setattr(yearlyout, "author", "Tammas Loughran")
    setattr(yearlyout, "contact", "t.loughran@student.unsw.edu.au")
    setattr(yearlyout, "date", dt.datetime.today().strftime('%Y-%m-%d'))
    setattr(yearlyout, "script", "ehfheatwaves.py")
    setattr(yearlyout, "dataset", "AWAP 0.5deg")
    setattr(yearlyout, "base_period", "%s-%s"%(str(bpstart),str(bpend)))
    setattr(yearlyout, "percentile", "%sth"%(str(pcntl)))
    otime = yearlyout.createVariable('time', 'f8', 'time', 
            fill_value=-999.99)
    setattr(otime, 'units', 'year')
    olat = yearlyout.createVariable('lat', 'f8', 'lat')
    setattr(olat, 'Longname', 'Latitude')
    setattr(olat, 'units', 'degrees_north')
    olon = yearlyout.createVariable('lon', 'f8', 'lon')
    setattr(olon, 'Longname', 'Longitude')
    setattr(olon, 'units', 'degrees_east')
    otpct = yearlyout.createVariable('t%spct'%(pcntl), 'f8', ('day','lat','lon'), 
            fill_value=-999.99)
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
    setattr(HWTout, 'Longname', 'First heat wave day of the year')
    otime[:] = range(syear, eyear)
    olat[:] = lats
    olon[:] = lons
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
    dailyout = Dataset('EHF_heatwaves_1911-2014_daily.nc', mode='w')
    dailyout.createDimension('time', len(time))
    dailyout.createDimension('lon', len(lons))
    dailyout.createDimension('lat', len(lats))
    setattr(dailyout, "author", "Tammas Loughran")
    setattr(dailyout, "contact", "t.loughran@student.unsw.edu.au")
    setattr(dailyout, "date", dt.datetime.today().strftime('%Y-%m-%d'))
    setattr(dailyout, "script", "ehfheatwaves.py")
    setattr(dailyout, "base_period", "%s-%s"%(str(bpstart),str(bpend)))
    setattr(dailyout, "percentile", "%sth"%(str(pcntl)))
    setattr(dailyout, "dataset", "AWAP 0.5deg")
    otime = dailyout.createVariable('time', 'f8', 'time',
                    fill_value=-999.99)
    setattr(otime, 'units', 'day as %Y%m%d.%f')
    setattr(otime, 'calendar', 'proleptic_gregorian')
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
    otime[:] = time
    olat[:] = lats
    olon[:] = lons
    dummy_array = np.ones((EHF.shape[0],)+original_shape[1:])*np.nan
    dummy_array[:,mask] = EHF
    oehf[:] = dummy_array.copy()
    dummy_array[:,mask] = event
    oevent[:] = dummy_array.copy()
    dummy_array[:,mask] = ends
    oends[:] = dummy_array.copy()
    dailyout.close()
