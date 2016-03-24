import sys
import numpy as np
import get_model_data as gmp
from netCDF4 import Dataset
from datetime import datetime, timedelta
import time
import utils

###########################################################################
#
#   AERIoe Model Observation Script
#   Author: Greg Blumberg  OU/CIMMS
#   Date Created: 1/28/2016
#
#   Recent changes to AERIoe by Dave Turner have included the option to 
#   have a thermodynamic profile be a "Y-vector" observation in the retrieval.
#   This facilitates the addition of radiosonde/Raman Lidar/DIAL data as inputs
#   into the AERIoe retrieval.  As a consequence, we have decided to remove the
#   AERI model prior concept from our processing and instead use a model thermodynamic
#   profile as another "observation" into the AERIoe retrieval.  This enables a better
#   retrieval in the mid-to-upper atmosphere.  This script will perform the conversions
#   from the model grids into the model-profile netCDF files that AERIoe can now
#   accept as inputs. 
#   
#   This script takes in a set of either real-time or archived model data
#   (e.g. ARM-RUC grids or realtime RAP/HRRR data or NCDC archived RUC data
#   or realtime GFS) and a lat/lon point.  It then pulls the temperature
#   and water vapor mixing ratio profile for each time at that point to develop 
#   a netCDF file.  The output file contains the time/height data for temperature and   
#   water vapor mixing ratio, pressure, height, temp_std, and wvmr_std.
#   It also contains metadata about how the file was made.
#
#   This code currently supports making these files from:
#   - realtime RAP profiles
#   - ARM-formatted RUC/RAP files
#
#   Arguments:
#       [1] YYYYMMDD - date of AERI file to run (ex: 20130531)
#       [2] VIP file - path to config file (ex: nwc_vip.txt)
#       [3] BHOUR - beginning hour of profiles (ex: 5)
#       [4] EHOUR - ending hour of profiles (ex: 5)
#
#   Copyright 2016
#
###########################################################################

def findVIPVariable(variable, filename):
    #This function searches for the value associatied with the key
    #(key is variable) within the VIP file and returns it either as a
    #float or a string
    try:
        config_file=open(filename,'r')
        ffile = config_file.read()
        config_file.close()
        ini=ffile.find(variable)+(len(variable)+1)
        rest=ffile[ini:]
        search_enter=rest.find('\n')
        data = (rest[:search_enter])
        datas = data.split(";")
        string = datas[0].replace('=', '')
    except:
        print "Could not read variable: \"" + variable + '\" from file: ' + filename
        sys.exit()

    try: 
        return float(string)
    except:
        return string
    
print "Starting generation of model grid \"observation\" files for use in AERIoe."

# Exit the program if the user didn't put in enough command line arguments
if len(sys.argv) != 5:
    print "Too few arguments in the command line, aborting."
    sys.exit() 
    
# Load in the command line arguments.
yyyymmdd = sys.argv[1] #Need date we're going to generate the observation.
vip = sys.argv[2] # Need the VIP file to figure out latitude/longitude point and what type of data source we're dealing with.
bhour = sys.argv[3] # What should the beginning hour of the output file be?
ehour = sys.argv[4] # What should the ending hour of the output file be?

# Exit the program if the YYYYMMDD argument isn't formatted correctly.
if len(yyyymmdd) != 8:
    print "YYYYMMDD argument not formatted correctly. Aborting."
    sys.exit()

 
print "Reading in the VIP variables..."
# Hardcoding these for debugging.
data_source = findVIPVariable('data_source', vip)
lat = findVIPVariable("aeri_lat", vip)
lon = findVIPVariable("aeri_lon", vip)
temporal_mesh_size = findVIPVariable('temporal_mesh_size', vip)
spatial_mesh_size = findVIPVariable('spatial_mesh_size', vip)

### Select the data source to generate the observation files.
if data_source == 1:
    print "VIP file says to generate model observation files using the 13 km RAP historical data from the NOAA NOMADS server."
    print "WARNING! Files are sometimes missing in random places on this server."
elif data_source == 2:
    print "VIP file says to generate model observation files using ARM-formatted RUC/RAP files."
elif data_source == 3:
    print "VIP file says to generate model observation files using the RAP MOTHERLODE UCAR datastream."
    print "This datastream contains more recent data and is recommended for realtime AERIoe runs."
else:
    print "Invalid value for \"data_source\" variable in the VIP file, aborting program."
    sys.exit()

# Convert the date/time strings into datetime objects for easier date/time manipulation
begin_dt = datetime.strptime(yyyymmdd + bhour, '%Y%m%d%H')
end_dt = datetime.strptime(yyyymmdd + ehour, '%Y%m%d%H')

print "Going to generate a model observation file for the times of: " + datetime.strftime(begin_dt, '%Y-%m-%d %H UTC') + ' to ' + datetime.strftime(end_dt, '%Y-%m-%d %H UTC')

if data_source == 1:
    # Using RUC/RAP historical data from the NOAA NOMADS server."
    output = gmp.getNOMADSModelObs(begin_dt, end_dt, temporal_mesh_size, spatial_mesh_size, lon, lat)
    output['arm_model_dir'] = 'n/a'
elif data_source == 2:
    # Most often used for the 1998-2003 ARM Boundary Facilities dataset
    # Using ARM-formatted RUC/RAP files to generate the model "observations"
    arm_model_dir = findVIPVariable('arm_model_dir', vip)
    output = gmp.getARMModelObs(arm_model_dir, begin_dt, end_dt, temporal_mesh_size, spatial_mesh_size, lon, lat) 
    output['arm_model_dir'] = arm_model_dir
elif data_source == 3:
    # Using the MOTHERLODE UCAR datasets when getting realtime observations.
    output = gmp.getRealtimeProfiles(begin_dt, end_dt, temporal_mesh_size, spatial_mesh_size, lon, lat)
    output['arm_model_dir'] = 'n/a'

output['spatial_mesh_size'] = spatial_mesh_size
output['temporal_mesh_size'] = temporal_mesh_size
output['aeri_lat'] = lat
output['aeri_lon'] = lon

# Get the path to write the converted model grid data into.
output['output_dir'] = findVIPVariable('output_dir', vip).strip()

# Performs a unit check on the data that is being save and converts it if it fails.
output = utils.unitCheck(output)

# Output the converted data into a netCDF file.
gmp.makeFile(output)

print "DONE."

