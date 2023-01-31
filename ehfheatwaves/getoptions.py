# -*- coding: utf-8 -*-
"""
getoptions.py parses the command line arguments to an options object.

Created on Fri Apr 13 21:58:40 2018

@author: Tammas Loughran
"""
import sys
import warnings
from optparse import OptionParser

options = None


class NoTmaxTminFileError(Exception):
    """Exception to be raised when no input data files are provided.
    """

    def __init__(self):
        print("Please specify tmax and tmin files.")


class InvalidBPFormatError(Exception):
    """Exception to be raised if the base period is not format is not ????-????
    """

    def __init__(self, baseperiod):
        print("The provided base period (", baseperiod, ")  must be formatted ????-????)")


class InvalidSeasonError(Exception):
    """Exception to be raised when the  given season is not summer or winter.
    """

    def __init__(self, season):
        print("The provided season much be either winter or summer. You specified ", season)


def parse_arguments(arguments):
    """parse_arguments parses the arguments to an options object, and handles some errors."""
    # Construct the options for the parser
    parser = OptionParser(usage="usage: %prog -x <FILE> -n <FILE> [options]")
    parser.add_option('-x', '--tmax', dest='tmaxfile', help='file containing tmax', metavar='FILE')
    parser.add_option('--vnamex', dest='tmaxvname', default='tasmax', help='tmax variable name', metavar='STR')
    parser.add_option('-n', '--tmin', dest='tminfile', help='file containing tmin', metavar='FILE')
    parser.add_option('--vnamen', dest='tminvname', default='tasmin', help='tmin variable name', metavar='STR')
    parser.add_option('--bpfx', dest='bpfx', help=('Indicates a future simulation, specifying a tmax file containing the historical base period to be used'), metavar='FILE')
    parser.add_option('--bpfn', dest='bpfn', help=('Indicates a future simulation, specifying a tmin file containing the historical base period to be used'), metavar='FILE')
    parser.add_option('-m', '--mask', dest='maskfile', help='file containing land-sea mask', metavar='FILE')
    parser.add_option('--vnamem', dest='maskvname', default='sftlf', help='mask variable name', metavar='STR')
    parser.add_option('--vnamet', dest='timevname', default='time', help='time variable name', metavar='STR')
    parser.add_option('-s', '--season', dest='season', default='summer', help='Season for annual metrics. Defaults to summer', metavar='STR')
    parser.add_option('-p', dest='pcntl', type='float', default=90, help='the percentile to use for thresholds. Defaults to 90', metavar='INT')
    parser.add_option('--base', dest='bp', default='1961-1990', help='base period to calculate thresholds. Default 1961-1990', metavar='YYYY-YYYY')
    parser.add_option('-q', '--qmethod', dest='qtilemethod', default='climpact', help='quantile interpolation method. Default is climpact', metavar='STR')
    parser.add_option('-d', '--daily', action="store_true", dest='daily', default=False, help='output daily EHF values and heatwave indicators')
    parser.add_option('--ehi', dest='ehi', action='store_true', default=False, help='Save the EHI values')
    parser.add_option('--dailyonly', action="store_true", dest='dailyonly', help='output only daily EHF values and suppress yearly output')
    parser.add_option('--t90pc', action="store_true", dest='t90pc', help='Calculate tx90pc and tn90pc heatwaves')
    parser.add_option('--tx90pc', action="store_true", dest='tx90pc', help='Calculate tx90pc seasonal heatwaves')
    parser.add_option('--tn90pc', action="store_true", dest='tn90pc', help='Calculate tn90pc seasonal heatwaves')
    parser.add_option('--tx90pc-daily', action="store_true", dest='tx90pcd', help='Calculate tx90pc daily heatwaves')
    parser.add_option('--tn90pc-daily', action="store_true", dest='tn90pcd', help='Calculate tn90pc daily heatwaves')
    parser.add_option('--noehf', action="store_true", dest='noehf', help='Supress EHF output and only use the specified t90pc')
    parser.add_option('-v', action="store_true", dest='verbose', help='Verbose')
    parser.add_option('--old-method', action="store_true", dest='oldmethod', help='Use the old definition of within-season heatwaves')
    parser.add_option('--invert-mask', action='store_true', dest='invertmask', help='Invert the land-sea mask.')
    parser.add_option('--flip-mask', action='store_true', dest='flipmask', help='Flip the mask upside down by latitude.')

    # Parse command line arguments to the options object
    global options
    options, args = parser.parse_args(arguments)

    # Handle errors and warnings
    if not options.tmaxfile and not options.tminfile: raise NoTmaxTminFileError
    if not '-' in options.bp: raise InvalidBPFormatError(options.bp)
    try:
        int(options.bp[:options.bp.index('-')])
        int(options.bp[options.bp.index('-')+1:])
    except ValueError:
        print('Base period years are not numbers.')
    assert int(options.bp[:4])<int(options.bp[5:9]), "Base period start is after end year."
    if (options.season!='summer')&(options.season!='winter'): raise InvalidSeasonError(options.season)
    warnmsg = "You didn't specify a land-sea mask. It's faster if you do, so this might take a while."
    if not options.maskfile: warnings.warn(warnmsg, UserWarning)

    # Additional options
    options.bpstart = int(options.bp[:4])
    options.bpend = int(options.bp[5:9])
    if options.daily or options.tx90pcd or options.tn90pcd or options.dailyonly:
        options.dailyout = True
        options.yearlyout = True
        if options.dailyonly: options.yearlyout = False
    else:
        options.dailyout = False
        options.yearlyout = True
    if options.t90pc:
        options.tx90pc = True
        options.tn90pc = True
    options.keeptmax = False
    if options.tx90pc or options.tx90pcd: options.keeptmax = True
    options.keeptmin = False
    if options.tn90pc or options.tn90pcd: options.keeptmin = True
    options.keeptave = True
    if options.noehf: options.keeptave = False

    return options
