# ehfheatwaves
A script to calculate heatwaves from AWAP and CMIP5 gridded datasets.

The original intention of this program was to provide daily EHF values and
heatwave indicators for any season rather than just annual or seasonal 
heatwave aspects. The script has been improved to handle CMIP5 quirks 
such different calendars.

Usage: ehfheatwaves_CMIP5.py -x <FILE> -n <FILE> -m <FILE> [options]

Options:
  -h, --help            show this help message and exit
  -x FILE, --tmax=FILE  file containing tmax
  --vnamex=STR          tmax variable name
  -n FILE, --tmin=FILE  file containing tmin
  --vnamen=STR          tmin variable name
  -m FILE, --mask=FILE  file containing land-sea mask
  --vnamem=STR          mask variable name
  -s STR, --season=STR  austal season for annual metrics. Defaults to austral
                        summer
  -p INT                the percentile to use for thresholds. Defaults to 90
  --base=YYYY-YYYY      base period to calculate thresholds. Default 1961-1990
  -q STR, --qmethod=STR
                        quantile interpolation method. Default is climpact
  -d, --daily           output daily EHF values and heatwave indicators
  --dailyonly           output only daily values and suppress yearly output

Required packages (prefereably the latest versions)
 * numpy
 * datetime
 * pandas
 * netCDF4
 * netcdftime
 * optparse

You must provide a tasmax file, tasmin file and a land sea mask. 
e.g. 
python ehfheatwaves_CMIP5.py -x "/direcory/tasmax_*" -n "/direcory/tasmin_*" -m "/direcory/mask.nc" -d
Make sure to use quotes if specifying multiple files with wild cards

Be careful when using the -d flag. This writes daily output of EHF values 
and heatwave indicators. During processing, unused data are removed using
the mask. When saving, this data need to be placed back into a very large
gridded array. This takes up a lot of ram ~30-40GB, so make sure that you
have enough ram.

To do list:
 * Save EHIs to file
 * Calculate heatwaves based on TX90pct
 * Save HWA using temp rather than the EHF
