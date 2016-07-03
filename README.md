# ehfheatwaves
A script to calculate heatwaves from AWAP and CMIP5 gridded datasets.

The original intention of this program was to provide daily EHF values and
heatwave indicators for any season rather than just annual or seasonal
heatwave aspects. The script has been improved to handle CMIP5 quirks
such different calendars.

Usage: ehfheatwaves.py -x <FILE> -n <FILE> -m <FILE> [options]

Options:
  -h, --help            show this help message and exit
  -x FILE, --tmax=FILE  file containing tmax
  --vnamex=STR          tmax variable name
  -n FILE, --tmin=FILE  file containing tmin
  --vnamen=STR          tmin variable name
  --bpfx=FILE           Indicates a future simulation, specifying a tmax file
                        containing the historical base period to be used
  --bpfn=FILE           Indicates a future simulation, specifying a tmin file
                        containing the historical base period to be used
  -m FILE, --mask=FILE  file containing land-sea mask
  --vnamem=STR          mask variable name
  --vnamet=STR          time variable name
  -s STR, --season=STR  Season for annual metrics. Defaults to summer
  -p INT                the percentile to use for thresholds. Defaults to 90
  --base=YYYY-YYYY      base period to calculate thresholds. Default 1961-1990
  -q STR, --qmethod=STR
                        quantile interpolation method. Default is climpact
  -d, --daily           output daily EHF values and heatwave indicators
  --dailyonly           output only daily EHF values and suppress yearly
                        output
  --t90pc               Calculate tx90pc and tn90pc heatwaves
  --tx90pc              Calculate tx90pc seasonal heatwaves
  --tn90pc              Calculate tn90pc seasonal heatwaves
  --tx90pc-daily        Calculate tx90pc daily heatwaves
  --tn90pc-daily        Calculate tn90pc daily heatwaves
  --noehf               Supress EHF output and only use the specified t90pc
  -v                    Verbose
  --old-method          Use the old definition of within-season heatwaves


Required packages (prefereably the latest versions)  
 * numpy
 * datetime
 * pandas
 * netCDF4
 * netcdftime
 * optparse

Notes:  
You must provide a tasmax file, tasmin file. A land sea mask is optional but recommended.  
e.g. python ehfheatwaves.py -x "/direcory/tasmax_\*" -n "/direcory/tasmin_*" -m "/direcory/mask.nc" -d  
Make sure to use quotes if specifying multiple files with wild cards

Be careful when using the -d flag. This writes daily output of EHF values 
and heatwave indicators to a sepatate file. During processing, unused data 
are removed using the mask. When saving, this data need to be placed back 
into a very large gridded array. This takes up a lot of ram ~30-40GB, so 
make sure that you have enough ram.

TX90pct and TN90pct heatwaves are also available. Daily values are available but only output one at a time. 
These will be output to separate files.

Global datasets are treated differently depending on which hemisphere the 
data are in. The data are split by hemispheres, and summer heatwaves are
calculated for the summer of that hemisphere. Heatwaves that are close to 
the equator should be cautiously considered. There is a discontinuity of 
season at the equator, with values at latitude 0 treated as part of the
southern hemisphere. Additionaly, there is little seasonal variation in 90th 
percentile in the tropics and is sometimes somewhat bimodal. This raises 
questions about the intepretation of tropical heatwaves using these definitions.
Bear this in mind.

The output file names are automatically generated for the time being. If your 
input dataset does not have any metadata about which model simulation it is, 
then the filename will contain underscores. For the time being just manually rename
the file to reflect the input dataset. Ill deal with this soon. 

This script is not parallelized.

To do list:
 * Save EHIs to file
 * Use standard metadata conventions
