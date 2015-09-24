import warnings
warnings.filterwarnings("ignore")
import sys
try: 
    import pandas as pd
except ImportError:
    print "Please install pandas"
    sys.exit()
import numpy as np
import datetime as dt
import qtiler
from netCDF4 import Dataset

# climpact/python/matlab
qtilemethod = 'climpact'
# season (winter/summer)
season = 'summer'
# save daily EHF output
dailyout = True
yearlyout = True
# base period
bpstart = 1961
bpend = 1990
# calendar
calendar = 'gregorian'
if calendar=='gregorian' or '365':
    daysinyear = 365
    if season=='winter':
        seasonlen = 153
        startday = 121
        endday = 274
    else:
        seasonlen = 151
        startday = 304
        endday = 455
elif calendar=='360':
    daysinyear = 360
    seasonlen = 150
    if season=='winter':
        startday = 121
        endday = 271
    else:
        startday = 301
        endday = 451
# percentile
pcntl = 90

# Load data
tmax_fname = ('/srv/ccrc/data35/z5032520/AWAP/daily/tmax/'
        'AWAP_TX_1911-2014_0.5deg.nc')
tmin_fname = ('/srv/ccrc/data35/z5032520/AWAP/daily/tmin/'
        'AWAP_TN_1911-2014_0.5deg.nc')
tmaxnc = Dataset(tmax_fname)
tmax = tmaxnc.variables['tmax'][:]
tmaxnc.close()
tminnc = Dataset(tmin_fname)
tmin = tminnc.variables['tmin'][:]
tave = (tmax + tmin)/2
# Mask the data
masknc = Dataset(
            '/srv/ccrc/data35/z5032520/AWAP/mask/AWAP_Land-Sea-Mask_0.5deg.nc')
mask = masknc.variables['LSM'][:]
mask = mask.astype(np.bool)
masknc.close()
original_shape = tave.shape
tave = tave[...,mask]
del tmin
del tmax
time = tminnc.variables['time'][:]
lats = tminnc.variables['lat'][:]
lons = tminnc.variables['lon'][:]
tminnc.close()

# Calculate date range
syear = int(str(time[0])[:4])
eyear = int(str(time[-1])[:4])
start = dt.datetime(syear,
        int(str(time[0])[4:6]),
        int(str(time[0])[6:8]))
end = dt.datetime(eyear,
        int(str(time[-1])[4:6]),
        int(str(time[-1])[6:8]))
if calendar=='360':
    class calendar360:
        def __init__(self,syear,eyear):
            self.year = np.repeat(range(syear,eyear+1), 360, 0)
            self.month = np.tile(np.repeat(range(1,12+1), 30, 0), 104)
            self.day = np.tile(np.tile(range(1,30+1), 12), 104)
    dates = calendar360(syear, eyear)
else:
    dates = pd.date_range(start,end)

# Cut out base period
tave = tave[(dates.month!=2)|(dates.day!=29),...]
time = time[(dates.month!=2)|(dates.day!=29)]
dates = dates[(dates.month!=2)|(dates.day!=29)]
base = ((bpstart<=dates.year)&(dates.year<=bpend))
tave_base = tave[base,...]

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

# Calculate EHF
EHF = np.ones(tave.shape)*np.nan
for i in xrange(tave.shape[0]):
    if i<31: continue
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

# Identify when heatwaves start with duration
# Given that first day contains duration
diff = np.zeros(event.shape)
diff[1:,...] = np.diff(event, axis=0)
ends = np.zeros(tave.shape,dtype=np.int)
ends[diff>2] = event[diff>2]

# Remove events less than 3 days
event[diff==2] = 0
event[np.roll(diff==2, 1, axis=0)] = 0
event[diff==1] = 0
del diff
event[event>0] = 1
event = event.astype(np.bool)
ends[ends<3] = 0

# Calculate metrics year by year
nyears = len(range(syear,eyear))
HWN = np.ones((nyears,tave.shape[1]))*np.nan
HWF = HWN.copy()
HWD = HWN.copy()
HWA = HWN.copy()
HWM = HWN.copy()
HWT = HWN.copy()
for iyear, year in enumerate(xrange(syear,eyear)):
    if (year==eyear)&(season=='summer'): continue 
    HWMtmp = np.ones((daysinyear,tave.shape[1]))*np.nan
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

if yearlyout:
    # Save to netCDF
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
