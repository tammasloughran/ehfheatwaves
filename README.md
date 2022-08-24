# ehfheatwaves
A tool to calculate heatwaves from gridded daily datasets.

[![DOI](https://zenodo.org/badge/41076341.svg)](https://zenodo.org/badge/latestdoi/41076341)

The default definition of a heatwave is any hot event of at least three days
duration where the Excess Heat Factor index is greater than 0.
One could alternatively use the tx90pc index (Tmax greater than the 90th
percentile) or tn90pc (Tmin greater than the 90th percentile.)

Seaonal heatwave aspect statistics are saved to a 'yearly' output netCDF file,
or more granular daily indicators can be saved to a 'daily' file.

The original intention of this program was to provide daily EHF values and
heatwave indicators for any season rather than just annual or seasonal
heatwave aspects. The script has since been expanded and improved to handle
CMIP5 quirks (such different calendars), and use alternative heatwave indices.

## Installation

Optionaly set up a virtualenv:
```
virtualenv heatwave_env
source heatwave_env/bin/activate
```

Install:
```
git clone https://github.com/tammasloughran/ehfheatwaves
cd ehfheatwaves
pip install .
```

## Usage

```
Usage: ehfheatwaves -x <FILE> -n <FILE> [options]

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
  --ehi                 Save the EHI values
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
```

## Notes

If you publish any research or other work that used this software, PLEASE cite
myself and this github repository in the acknowledgments section of your
article. Send the final print to me as I would also be interested in reading
your work. Thank you.


Input and output files:

You must provide at least some temperature data using the -x or -n arguments.
If using EHF, you must provide both tmax and tmin (-x and -n). If only using
tx90pc or tn90pc indices then you may only provide the relevant data for each.
Specifying a land sea mask is optional, using the -m argument, but recommended.
Multiple files can be specified using wildcards, use quotes when doing so. e.g.
```
python ehfheatwaves.py -x "/direcory/tasmax_\*" -n "/direcory/tasmin_*" -m "/direcory/mask.nc" -d
```

Be careful when using the -d flag. This writes daily output of EHF values
and heatwave indicators to a sepatate file. During processing, unused data
are removed using the mask. When saving, this data need to be placed back
into a very large gridded array. This might take up a lot of ram ~30-40GB, so
make sure that you have enough if running on a desktop.

The output file names are automatically generated. Output files are named:
```
<definition>_heatwaves_<modelname>_<experiment>_<ripcode>_<frequency>.nc
```
definition = EHF/tx90pct/tn90pct

modelname = String taken from the 'model_id' attribute in your input netCDF file.

experiment = String taken from the 'parent_experiment_rip' attribute.

ripcode = rXiXpX taken from 'realization', 'initialization_method' and ''physics_version' attributes.

frequency = yearly/daily

If your input dataset does not have any metadata about which model simulation
it is, then the filename will just contain underscores. If you are running the
program several times then make sure to rename the output files so they are not
overwritten.

More info on output data:

Global datasets are treated differently depending on which hemisphere the
data are in. The data are split by hemispheres, and summer heatwaves are
calculated for the summer of that hemisphere. Heatwaves that are close to
the equator should be treated with caution. There is a discontinuity of
season at the equator, with values at latitude 0 treated as part of the
southern hemisphere. Additionaly, there is little seasonal variation in 90th
percentile in the tropics and is sometimes somewhat bimodal. This raises
questions about the intepretation of tropical heatwaves using these definitions.
Bear this in mind.

There are three daily indocators. 'EHF'/'tx90pct'/'tn90pct' are the raw index
values for where heatwaves occur. They are 0 if a heatwave is not occuring.
'events' is a simple boolean indicator for whether a heatwave is occuring on
any given day.
'ends' denotes the duration of a heatwave on the first day. All other heatwave
days are 0.

To do list:
 * Use standard metadata conventions

