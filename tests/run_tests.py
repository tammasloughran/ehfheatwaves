#!/usr/bin python3
# -*- coding: utf-8 -*-
"""
run_tests.py tests the accuracy of ehfheatwaves.py compared to heatwave data
calculated by climpact2 R package, and runs some basic unit tests.

@author: Tammas Loughran
"""

# Load modules
import os
import netCDF4 as nc
import sys
sys.path.append('../')
import numpy as np
import unittest
from ehfheatwaves import *
import qtiler


class TestRQtiler(unittest.TestCase):
    """Test the R based quantile function"""

    testdata = np.array([1,2,3,4,5,6,7,8,9])

    def testNegativeP(self):
        self.assertRaises(ValueError,qtiler.quantile_R(self.testdata,-1))

    def testZeroP(self):
        self.assertEqual(qtiler.quantile_R(self.testdata,0), min(self.testdata))

    def testHundredP(self):
        self.assertEqual(qtiler.quantile_R(self.testdata,100), max(self.testdata))

    def testHundredPlusP(self):
        self.assertRaises(ValueError,qtiler.quantile_R(self.testdata,101))

    def testFractionOneP(self):
        self.assertEqual(qtiler.quantile_R(self.testdata,1,fraction=True), max(self.testdata))

    def testPercentAsFraction(self):
        self.assertRaises(ValueError,qtiler.quantile_R(self.testdata,50,fraction=True))

    def testFractionAsPercent(self):
        self.assertRaises(Warning,qtiler.quantile_R(self.testdata,0.5,fraction=False))

    def testInvalidItype(self):
        self.assertRaises(ValueError,qtiler.quantile_R(self.testdata,50,itype=0))

    def testKnownCases(self):
        knownCases = ((1,7),(2,7),(3,6),(4,6.3),(5,6.8),(6,7),(7, 6.6), (8, 6.866667), (9,6.85))
        for case in knownCases:
            self.assertAlmostEqual(qtiler.quantile_R(self.testdata,70,itype=case[0]),case[1],places=6)

class TestZhangQtiler(unittest.TestCase):
    """Test the Zhang quantile function"""

    testdata = np.array([1,2,3,4,5,6,7,8,9])

    def testNegativeP(self):
        self.assertRaises(ValueError,qtiler.quantile_zhang(self.testdata,-1))

    def testZeroP(self):
        self.assertEqual(qtiler.quantile_zhang(self.testdata,0), min(self.testdata))

    def testHundredP(self):
        self.assertEqual(qtiler.quantile_zhang(self.testdata,100), max(self.testdata))

    def testHundredPlusP(self):
        self.assertRaises(ValueError,qtiler.quantile_zhang(self.testdata,101))

    def testFractionOneP(self):
        self.assertEqual(qtiler.quantile_zhang(self.testdata,1,fraction=True), max(self.testdata))

    def testPercentAsFraction(self):
        self.assertRaises(ValueError,qtiler.quantile_zhang(self.testdata,50,fraction=True))

    def testFractionAsPercent(self):
        self.assertRaises(Warning,qtiler.quantile_zhang(self.testdata,0.5,fraction=False))

    def testInvalidItype(self):
        self.assertRaises(ValueError,qtiler.quantile_zhang(self.testdata,50,itype=0))

    def testKnownCase(self):
        self.assertEqual(qtiler.quantile_zhang(self.testdata,70),7.0)


if __name__=='__main__':
    # First test the script as a whole against climpact2 data.
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

    unittest.main()
