"""run_tests.py tests the accuracy of ehfheatwaves.py compared to heatwave data
calculated by climpact2 R package.
"""
# Load modules
import os
import netCDF4 as nc
import sys
import numpy as np

# Define command and arguments
run_comnd = 'python ehfheatwaves.py'
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
    print error_msg
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
difference = thwdata - chwdata

# Compare difference to a threshold
threshold = 0.00001
if np.any(difference>threshold):
    print "Fail"
else:
    print "Pass"

# Remove the test file
os.remove('testfile.nc')