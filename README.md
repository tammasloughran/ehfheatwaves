# ehfheatwaves
A script to calculate heatwaves from AWAP and CMIP5 gridded datasets.

The original intention of this program was to provide daily EHF values and
heatwave indicators for any season rather than just annual or seasonal 
heatwave aspects. The script has been improved to handle CMIP5 quirks 
such different calendars.

You must provide a tasmax file, tasmin file and a land sea mask. 
e.g. 
python ehfheatwaves_CMIP5.py -x '/direcory/tasmax_*' -n '/direcory/tasmin_*' -m '/direcory/mask.nc' -d

Be careful when using the -d flag. This writes daily output of EHF values 
and heatwave indicators. During processing, unused data are removed using
the mask. When saving, this data need to be placed back into a very large
gridded array. This takes up a lot of ram ~70GB, so make sure that you 
are the only user on your server.
