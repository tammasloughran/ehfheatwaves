#!/usr/bin python3
# -*- coding: utf-8 -*-
"""
run_tests.py tests the accuracy of ehfheatwaves.py compared to heatwave data
calculated by climpact2 R package, and runs some unit tests.

@author: Tammas Loughran
"""

# Load modules
import os
import netCDF4 as nc
import sys
sys.path.append('../')
import numpy as np
import datetime as dt
import unittest
from ehfheatwaves import *
from getoptions import *
import optparse
import qtiler


class TestRQtiler(unittest.TestCase):
    """Test the R based quantile function."""

    testdata = np.array([1,2,3,4,5,6,7,8,9])

    def testNegativeP(self):
        """p < 0 should return an exception error."""
        self.assertRaises(ValueError, qtiler.quantile_R, self.testdata, -1)

    def testZeroP(self):
        """p = 0 should return the lowest value."""
        self.assertEqual(qtiler.quantile_R(self.testdata,0), min(self.testdata))

    def testHundredP(self):
        """p = 100 should return the largest value."""
        self.assertEqual(qtiler.quantile_R(self.testdata,100), max(self.testdata))

    def testHundredPlusP(self):
        """p > 100 should return an exception error."""
        self.assertRaises(ValueError, qtiler.quantile_R, self.testdata, 101)

    def testFractionOneP(self):
        """p = 1 as a fraction should return the largest value."""
        self.assertEqual(qtiler.quantile_R(self.testdata,1,fraction=True), max(self.testdata))

    def testPercentAsFraction(self):
        """p (when fraction) > 1 should return an exception Error."""
        self.assertRaises(ValueError, qtiler.quantile_R, self.testdata, 50, fraction=True)

    def testFractionAsPercent(self):
        """p (when percent) < 1 should give a warning."""
        self.assertWarns(Warning, qtiler.quantile_R, self.testdata, 0.5, fraction=False)

    def testInvalidItype(self):
        self.assertRaises(KeyError, qtiler.quantile_R, self.testdata, 50, itype=0)

    def testKnownCases(self):
        """Test known cases for corectness."""
        knownCases = ((1,7),(2,7),(3,6),(4,6.3),(5,6.8),(6,7),(7, 6.6), (8, 6.866667), (9,6.85))
        for case in knownCases:
            self.assertAlmostEqual(qtiler.quantile_R(self.testdata,70,itype=case[0]),case[1],places=6)


class TestZhangQtiler(unittest.TestCase):
    """Test the Zhang quantile function."""

    testdata = np.array([1,2,3,4,5,6,7,8,9])

    def testNegativeP(self):
        """p < 0 should return an exception error."""
        self.assertRaises(ValueError, qtiler.quantile_zhang, self.testdata, -1)

    def testZeroP(self):
        """p = 0 should return the lowest value."""
        self.assertEqual(qtiler.quantile_zhang(self.testdata,0), min(self.testdata))

    def testHundredP(self):
        """p = 100 should return the largest value."""
        self.assertEqual(qtiler.quantile_zhang(self.testdata,100), max(self.testdata))

    def testHundredPlusP(self):
        """p > 100 should return an exception error."""
        self.assertRaises(ValueError, qtiler.quantile_zhang, self.testdata, 101)

    def testFractionOneP(self):
        """p = 1 as a fraction should return the largest value."""
        self.assertEqual(qtiler.quantile_zhang(self.testdata,1,fraction=True), max(self.testdata))

    def testPercentAsFraction(self):
        """p (when fraction) > 1 should return an exception Error."""
        self.assertRaises(ValueError, qtiler.quantile_zhang, self.testdata, 50, fraction=True)

    def testFractionAsPercent(self):
        """p (when percent) < 1 should give a warning."""
        self.assertWarns(Warning, qtiler.quantile_zhang, self.testdata,  0.5, fraction=False)

    def testKnownCase(self):
        """Test a known case for corectness."""
        self.assertEqual(qtiler.quantile_zhang(self.testdata,70),7.0)


class TestClimpactQtiler(unittest.TestCase):
    """Test the climapact quantile function"""

    testdata = np.array([1,2,3,4,5,6,7,8,9])

    def testNegativeP(self):
        """p < 0 should return an exception error."""
        self.assertRaises(ValueError, qtiler.quantile_climpact, self.testdata, -1)

    def testZeroP(self):
        """p = 0 should return the lowest value."""
        self.assertEqual(qtiler.quantile_climpact(self.testdata,0), min(self.testdata))

    def testHundredP(self):
        """p = 100 should return the largest value."""
        self.assertEqual(qtiler.quantile_climpact(self.testdata,100), max(self.testdata))

    def testHundredPlusP(self):
        """p > 100 should return an exception error."""
        self.assertRaises(ValueError, qtiler.quantile_climpact, self.testdata, 101)

    def testFractionOneP(self):
        """p = 1 as a fraction should return the largest value."""
        self.assertEqual(qtiler.quantile_climpact(self.testdata,1,fraction=True), max(self.testdata))

    def testPercentAsFraction(self):
        """p (when fraction) > 1 should return an exception Error."""
        self.assertRaises(ValueError, qtiler.quantile_climpact, self.testdata, 50, fraction=True)

    def testFractionAsPercent(self):
        """p (when percent) < 1 should give a warning."""
        self.assertRaises(Warning, qtiler.quantile_climpact, self.testdata, 0.5, fraction=False)

    def testKnownCase(self):
        """Test a known case for corectness."""
        self.assertAlmostEqual(qtiler.quantile_climpact(self.testdata,70),6.8666666666666654)

class TestCalendar360(unittest.TestCase):
    """Tests for the Calendar360 class."""

    startdate = dt.datetime(1970,1,1)
    enddate = dt.datetime(2000,1,1)
    calendar = calendar360(startdate, enddate)

    def testEndBeforeStart(self):
        """Should return an exception if the end date is before the start date."""
        self.assertRaises(DatesOrderError, calendar360, self.enddate, self.startdate)

    def testAtributes(self):
        """Calendar object should contain 3 attributes contaning arrays."""
        self.assertIs(type(self.calendar.day), np.ndarray)
        self.assertIs(type(self.calendar.month), np.ndarray)
        self.assertIs(type(self.calendar.year), np.ndarray)

    def testIndexable(self):
        """Calendar object should be indexable, which returns a date."""
        self.assertIs(type(self.calendar[0]), dt.datetime)


class TestIdentifyHW(unittest.TestCase):
    """Tests for the identify_hw function."""

    # Some example EHF index data
    ehfdata = np.array([-1.3, -0.3, 3.4, 8.5, -0.4, # Not a heatwave
                        5.6, 10.2, 20.4, -1.4, # a three day heatwave starting on day 6
                        7.8, 15.5, 16.9, 17.9, 30.2, -3.3]) # a 5 day heatwave starting on day 10

    known_events = np.array([0,0,0,0,0,1,1,1,0,1,1,1,1,1,0])
    known_ends = np.array([0,0,0,0,0,3,0,0,0,5,0,0,0,0,0])

    def testReturnTupple(self):
        """Should return a tupple containging the event indicator and the durations (numpy.ndarray)."""
        result = identify_hw(self.ehfdata)
        self.assertIs(tuple, type(result))
        self.assertIs(np.ndarray, type(result[0]))
        self.assertIs(np.ndarray, type(result[1]))

    def testKnownEventsEnds(self):
        """Tests the known cases."""
        events, ends = identify_hw(self.ehfdata)
        self.assertTrue((events==self.known_events).all())
        self.assertTrue((ends==self.known_ends).all())

    def testShape(self):
        """The shape of the input EHF index and the outputs should be the same."""
        events, ends = identify_hw(self.ehfdata)
        input_shape = self.ehfdata.shape
        self.assertEqual(events.shape, input_shape)
        self.assertEqual(ends.shape, input_shape)


class TestIdentifySemiHW(unittest.TestCase):
    """Tests for the identify_semi_hw function."""

    # Some example EHF index data
    ehfdata = np.array([1.3, 0.3, -3.4, 8.5, -0.4, # 2 day heatwave
                        5.6, 10.2, 20.4, -1.4, # a three day heatwave starting on day 6
                        7.8, 15.5, 16.9, 17.9, 30.2, -3.3]) # a 5 day heatwave starting on day 10

    known_events = np.array([1,1,0,0,0,1,1,1,0,1,1,1,1,1,0])
    known_ends = np.array([2,0,0,0,0,3,0,0,0,5,0,0,0,0,0])

    def testReturnTupple(self):
        """Should return a tupple containging the event indicator and the durations (numpy.ndarrays)."""
        result = identify_hw(self.ehfdata)
        self.assertIs(tuple, type(result))
        self.assertIs(np.ndarray, type(result[0]))
        self.assertIs(np.ndarray, type(result[1]))

    def testKnownEventsEnds(self):
        """Tests the known cases."""
        events, ends = identify_hw(self.ehfdata)
        self.assertTrue((events==self.known_events).all())
        self.assertTrue((ends==self.known_ends).all())

    def testShape(self):
        """The shape of the input EHF index and the outputs should be the same."""
        events, ends = identify_hw(self.ehfdata)
        input_shape = self.ehfdata.shape
        self.assertEqual(events.shape, input_shape)
        self.assertEqual(ends.shape, input_shape)


class TestGetOptions(unittest.TestCase):
    """Tests for the getoptins module"""

    args = ['-x', 'tests/climpact2.sampledata.gridded.1991-2010.nc',
            '-n', 'tests/climpact2.sampledata.gridded.1991-2010.nc',
            '--vnamex=tmax' '--vnamen=tmin', '--base=1991-2010', '-v']

    def testConstructsOptions(self):
        """The parse_arguments function should return an options object"""
        self.assertIs(type(getoptions.parse_arguments(self.args)), optparse.Values)

    def testNoFilesError(self):
        """Should return an exception if no data files are provided."""
        self.assertRaises(NoTmaxTminFileError, getoptions.parse_arguments, ['-v'])

    def testInvalidBP(self):
        """Should return an exception if the base period is not in the correct format."""
        args = ['-x', 'tests/climpact2.sampledata.gridded.1991-2010.nc',
            '-n', 'tests/climpact2.sampledata.gridded.1991-2010.nc',
            '--base=19912010']
        self.assertRaises(InvalidBPFormatError, getoptions.parse_arguments, args)

    def testInvalidSeason(self):
        """Should return an exception if an invalis season is provided."""
        args = ['-x', 'tests/climpact2.sampledata.gridded.1991-2010.nc',
            '-n', 'tests/climpact2.sampledata.gridded.1991-2010.nc',
            '--base=1991-2010', '--season=autmn']
        self.assertRaises(InvalidSeason, getoptions.parse_arguments, args)

    def testBPOrder(self):
        """Should return an assertion error if the start year of the base period is before the end year"""
        args = ['-x', 'tests/climpact2.sampledata.gridded.1991-2010.nc',
            '-n', 'tests/climpact2.sampledata.gridded.1991-2010.nc',
            '--base=2010-1960']
        self.assertRaises(AssertionError, getoptions.parse_arguments, args)

    def testBPAreNums(self):
        """Should return a ValueError exception if the base period provided are not numbers"""
        args = ['-x', 'tests/climpact2.sampledata.gridded.1991-2010.nc',
            '-n', 'tests/climpact2.sampledata.gridded.1991-2010.nc',
            '--base=abc-def']
        self.assertRaises(ValueError, getoptions.parse_arguments, args)

    def testMaskWarning(self):
        """Should throw a warning if no land-sea mask is provided"""
        args = ['-x', 'tests/climpact2.sampledata.gridded.1991-2010.nc',
            '-n', 'tests/climpact2.sampledata.gridded.1991-2010.nc',
            '--base=abc-def']
        self.assertWarns(UserWarning, getoptions.parse_arguments, self.args)

if __name__=='__main__':
    # Test the script as a whole against climpact2 data.
    # Define command and arguments
    if sys.version[0]=='3':
        run_comnd = 'python3 ../ehfheatwaves.py'
    else:
        run_comnd = 'python ../ehfheatwaves.py'
    arguments = ' -x climpact2.sampledata.gridded.1991-2010.nc' \
        ' -n climpact2.sampledata.gridded.1991-2010.nc' \
        ' --vnamex=tmax' \
        ' --vnamen=tmin' \
        ' --base=1991-2010' \
        ' -v'

    # Check if outputfile already exists.
    file_exists = os.path.isfile('EHF_heatwaves____yearly_summer.nc')
    error_msg = 'The output files already exist. Please clear the output files '\
        'before running the test script.'
    if file_exists:
        print(error_msg)
        sys.exit(1)

    # Run ehfheatwaves.py
    os.system(run_comnd + arguments)

    # Rename the output file.
    os.rename('EHF_heatwaves____yearly_summer.nc', 'testfile.nc')

    # Load the results.
    testnc = nc.Dataset('testfile.nc','r')
    hwf = testnc.variables['HWF_EHF'][:]
    thwdata = np.ones((6,)+hwf.shape)
    thwdata[0,...] = hwf
    thwdata[1,...] = testnc.variables['HWN_EHF'][:]
    thwdata[2,...] = testnc.variables['HWD_EHF'][:]
    thwdata[3,...] = testnc.variables['HWA_EHF'][:]
    thwdata[4,...] = testnc.variables['HWM_EHF'][:]
    thwdata[5,...] = testnc.variables['HWT_EHF'][:]
    testnc.close()

    # Load the comparison data
    chwdata = np.load('hwdata.npy')

    # Take the difference between the test and comparison data
    difference = abs(thwdata - chwdata)

    # Compare difference to a threshold
    threshold = 0.00001
    if np.any(difference>threshold):
        print("Results accuracy: Fail")
    else:
        print("Results Accuracy: Pass")

    # Remove the test file
    os.remove('testfile.nc')

    # Execute unit tests
    unittest.main()