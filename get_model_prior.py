from netCDF4 import Dataset
#from pylab import *
import numpy as np
import sys
#from mpl_toolkits.basemap import Basemap
from datetime import datetime, timedelta
import glob
import platform
import utils

"""
    Script Name: get_model_prior.py
    
    Contains two primary parsing functions:
    getARMProfiles() - which will parse out the profiles within a specified spatial domain
                       for each ARM-formatted RUC/RAP model file.
    getMotherlodeProfiles() - will parse out the profiles contained within the 13 km RAP grid
                              (used for real-time model observations for AERIoe.)
                              
    Also contains two functions that call the parsing functions and control the files that 
    get opened in order to create the model observation files for AERIoe.
"""

def getARMProfiles(rap_path, yyyymmdd, hh,  aeri_lon, aeri_lat, size):
    """
    This function is used to parse through ARM-formatted RAP/RUC analysis files that can be downloaded
    from the ARM Archive.  These files are primarily used to provide improved upper-air information
    to the AERIoe retrieval, but for only AERI deployments that took place in the SGP.
    
    The files it looks for have the regular expression: *syn*yyyymmdd.hh*.cdf
                                                        *all*yyyymmdd.hh*.cdf
                                                        
    Using the yyyymmddhh string, this function searches for ARM-formatted RAP/RUC netCDF files.  
    Using the aeri_lat, aeri_lon this code will parse out profiles from within the user-specified
    spatial domain. 
    
    This parsing creates four arrays: temps, mxrs, press, and hghts
    These arrays are 2-D in the sense that the first index is the profile and the second index is 
    the vertical grid index of the thermodynamic profile.
    
    The parsing also returns the thermodynamic profile from the grid point closest to the AERI location.
    
    TODO: - Dump out the wind profile (u,v) information from this grid point too.
          - Ensure that this function works for different RAP/RUC resolution and files (e.g 60 km, 40 km)
    
    This function in the past has gotten called by getARMModelPrior(), which uses this function to open up
    several ARM-formatted RUC/RAP files across a certain time window.    
    
    Returns
    -------

    """   
    print "\tSearching for \"syn\" files."
    
    # Find files that have the "syn" in their filename
    files = glob.glob(rap_path.strip() + '/*syn*' + yyyymmdd + '.' + hh + '*.cdf') 
    
    # If you can't find any files with "syn" in their filename, search for those with the "all"
    if len(files) == 0:
        print "\tNone, searching for \"all\" files."
        files = glob.glob(rap_path.strip() + '/*all*' + yyyymmdd + '.' + hh + '*.cdf') 
    
    # Ensure that these variables are floats and ints since they may have been passed as strings.
    aeri_lon = float(aeri_lon)
    aeri_lat = float(aeri_lat)     
    size = int(size)

    # Ensure that the file found can be opened up.
    try:
         # Tell the user that data has been found
        print "\t\tFOUND DATA FOR PRIOR: " + files[0]
        d = Dataset(files[0])
    except:
        # Give messages to the user that the file couldn't be found or opened.
        # Return error codes.
        print "\t\tUnable to find ARM data to generate prior."
        print "\t\tARM RAP/RUC \"syn\" data needs to be in the directory: " + rap_path
        print "\t\tFor the date and hour of: " + yyyymmdd + ' ' + hh + ' UTC'
        return
        
    # Grab the latitude/longitude grid stored in the ARM-formatted model grid file
    lon = d.variables['longitude'][:]
    lat = d.variables['latitude'][:]
    
    # Find the indices associated with the grid point nearest to the instrument.
    idy, idx = utils.find_index_of_nearest_xy(lon, lat, aeri_lon, aeri_lat)
   
    # Extract the pressure, temperature, RH, and height grids for the specificed spatial domain
    # as well as the surface properties of these variables.  Profiles are on an isobaric grid.
    # Variables that begin with "center_" are for the profile nearest the AERI instrument.
    pres = d.variables['pressurepgrid'][:]
    temp =  d.variables['tempp'][0,idy-size:idy+size,idx-size:idx+size,:]
    center_temp =  d.variables['tempp'][0,idy,idx,:]
    rh = d.variables['rhp'][0,idy-size:idy+size,idx-size:idx+size,:]
    center_rh = d.variables['rhp'][0,idy,idx,:]
    hght = d.variables['heightgpp'][0,idy-size:idy+size,idx-size:idx+size,:]
    center_hght = d.variables['heightgpp'][0,idy,idx]
    
    sfc_pres = d.variables['pressuresrf'][0,idy-size:idy+size,idx-size:idx+size]/100. # Convert from Pa to mb
    center_sfc_pres = d.variables['pressuresrf'][0,idy,idx]/100. # Convert from Pa to mb
    sfc_hght = d.variables['heightsrf'][idy-size:idy+size,idx-size:idx+size]
    center_sfc_hght = d.variables['heightsrf'][idy,idx]
    sfc_temp = d.variables['temp2m'][0,idy-size:idy+size,idx-size:idx+size]
    center_sfc_temp = d.variables['temp2m'][0,idy,idx]
    sfc_rh = d.variables['rh2m'][0,idy-size:idy+size,idx-size:idx+size]
    center_sfc_rh = d.variables['rh2m'][0,idy,idx]
   
    # Find the indicies that correspond to the elements of the vertical grid that are above ground.
    idx_aboveground = np.where(center_sfc_hght < center_hght)[0]
        
    # Merge the 2 meter AGL variables with the rest of the above ground profile.
    center_hght = (np.hstack((center_sfc_hght+2, center_hght[idx_aboveground])) - center_sfc_hght)/1000.
    center_temp = np.hstack((center_sfc_temp, center_temp[idx_aboveground]))
    center_rh = np.hstack((center_sfc_rh, center_rh[idx_aboveground]))
    center_pres = np.hstack((center_sfc_pres, pres[idx_aboveground]))
    
    center_q = 1000. * utils.rh2q(center_temp, center_pres*100., center_rh/100.)
    
    # Close the ARM-formatted netCDF RAP/RUC file.
    d.close()
    
    # Initalize the profile storage arrays.
    mxrs = []
    temps = []
    press = []
    hghts = []
    
    # Loop over the horizontal grid and pull out the vertical profiles at each grid point.
    for index, x in np.ndenumerate(sfc_hght):
        # Find the indicies that correspond to the elements of the vertical grid that are above ground.
        idx_aboveground = np.where(sfc_hght[index] < hght[index[0], index[1],:])[0]
        
        # Merge the 2 meter AGL variables with the rest of the above ground profile.
        new_hght = (np.hstack((sfc_hght[index]+2, hght[index[0], index[1], idx_aboveground])) - sfc_hght[index])/1000.
        new_temp = np.hstack((sfc_temp[index], temp[index[0], index[1], idx_aboveground]))
        new_rh = np.hstack((sfc_rh[index], rh[index[0], index[1], idx_aboveground]))
        new_pres = np.hstack((sfc_pres[index], pres[idx_aboveground]))
        new_q = utils.rh2q(new_temp, new_pres*100., new_rh/100.)*1000.
        
        # Append the full ground to top of model profile to storage arrays.
        mxrs.append(new_q)
        temps.append(new_temp)
        press.append(new_pres)
        hghts.append(new_hght)
        
    # Convert the storage arrays to Numpy arrays
    temps = np.asarray(temps)
    mxrs = np.asarray(mxrs)
    press = np.asarray(press)
    hghts = np.asarray(hghts)
    
    # TODO: figure out the format of the data that gets returned.
    stop
    return

def getMotherlodeProfiles(yyyymmddhh, time_window, aeri_lat, aeri_lon, size=5):
    """
    This function is used to parse through RAP data located on the THREDDS Motherlode server
    
    Using the yyyymmddhh string, this function reads in the RAP data for that
    time from the Motherlode Server.  Using the aeri_lat, aeri_lon grids and an
    netCDF file of the surface terrain height @ 13 km and the time_window variable,
    this code will parse out RAP profiles from within this spatial and temporal mesh/data cube.
    This data contains both the analysis (at index=0) and the forecast data.  Analysis data is only
    included if time_window = 0.  If forecast data is to be used, time_window needs to be greater than 1.
    
    This parsing creates four arrays: temps, mxrs, press, and hghts
    These arrays are 2-D in the sense that the first index is the profile and the second index is 
    the vertical grid index of the thermodynamic profile.
    
    The parsing also returns the thermodynamic profile from the grid point closest to the AERI location.
    
    This function and methodology works best when we want to run AERIoe in realtime mode, as
    the RAP forecast and analysis data has a 2-week lifetime on the Motherlode server.
    
    TODO: - Include Unidata Siphon compatability (might work better!)
          - Dump out the wind profile (u,v) information from this grid point too.
    
    Returns
    -------

"""
    recent_rap_path = 'http://thredds.ucar.edu/thredds/dodsC/grib/NCEP/RAP/CONUS_13km/RR_CONUS_13km_' + \
        yyyymmddhh[:8] + '_' + yyyymmddhh[8:10] + '00.grib2/GC'
    print recent_rap_path
    try:
        d = Dataset(recent_rap_path)
        print "Found 13 km RAP data for this date on the Motherlode UCAR server."
        model_name = "RAP13km"
        path = recent_rap_path
    except:
        print "No RAP data found for this date on the Motherlode UCAR server."
        print "Data Path:", recent_rap_path
        return None
    
    # Read in the terrain data and the lat,lon grid for the 13 km RAP grid.
    ll = Dataset('13km_latlon.nc')
    lon = ll.variables['lon'][:]
    lat = ll.variables['lat'][:]
    
    # Ensure that the aeri_lat, aeri_lon variables are floats.
    aeri_lon = float(aeri_lon)
    aeri_lat = float(aeri_lat)
    
    # Find the indices for the nearest grid point to the AERI location
    idy, idx = utils.find_index_of_nearest_xy(lon, lat, aeri_lon, aeri_lat)

    # Read in the 13 km RAP data for parsing.
    pres = d.variables['isobaric'][:] #Pascals
    temp =  d.variables['Temperature_isobaric'][:time_window,:,idy-size:idy+size,idx-size:idx+size] # 
    center_temp =  d.variables['Temperature_isobaric'][:time_window,:,idy, idx] # 
    rh = d.variables['Relative_humidity_isobaric'][:time_window,:,idy-size:idy+size,idx-size:idx+size]
    center_rh = d.variables['Relative_humidity_isobaric'][:time_window,:,idy, idx]
    hght = d.variables['Geopotential_height_isobaric'][:time_window,:,idy-size:idy+size,idx-size:idx+size]
    center_hght = d.variables['Geopotential_height_isobaric'][:time_window,:,idy ,idx]
    sfc_pres = d.variables['Pressure_surface'][:time_window, idy-size:idy+size,idx-size:idx+size]  #Pascals
    center_sfc_pres = d.variables['Pressure_surface'][:time_window, idy, idx]  #Pascals
    sfc_hght = ll.variables['Geopotential_height_surface'][idx-size:idx+size, idy-size:idy+size]
    center_sfc_hght = ll.variables['Geopotential_height_surface'][idx, idy]
    #sfc_hght = d.variables['Geopotential_height_surface'][:time_window, idy-size:idy+size,idx-size:idx+size]
    sfc_temp = d.variables['Temperature_height_above_ground'][:time_window,0, idy-size:idy+size,idx-size:idx+size]
    center_sfc_temp = d.variables['Temperature_height_above_ground'][:time_window,0, idy, idx]
    sfc_rh = d.variables['Relative_humidity_height_above_ground'][:time_window,0, idy-size:idy+size,idx-size:idx+size]
    center_sfc_rh = d.variables['Relative_humidity_height_above_ground'][:time_window,0, idy, idx]
   
    # Merge the model profile data together to get the thermodynamic profile nearest to the AERI location.
    hght_prof = center_hght[0, :]
    idx_aboveground = np.where(pres < center_sfc_pres[0])[0]
    point_hght = np.hstack((center_sfc_hght+2, center_hght[0, idx_aboveground][::-1]))
    new_temp = np.hstack((center_sfc_temp[0], center_temp[0, idx_aboveground][::-1]))
    new_rh = np.hstack((center_sfc_rh[0], center_rh[0, idx_aboveground][::-1]))
    new_pres = np.hstack((center_sfc_pres[0], pres[idx_aboveground][::-1]))
    new_q = utils.rh2q(new_temp, new_pres, new_rh/100.)*1000.
   
    # Close the terrain dataset and the 13 km RAP dataset.
    ll.close()
    d.close()
    
    # Intitalize the profile storage arrays.
    mxrs = []
    temps = []
    press = []
    hghts = []
    
    # Loop over the time index for the RAP data.
    for t in range(len(temp)):
        # Loop over each point in the horizontal grid, merge the ground properties with the upper air data and store the profile.
        for index, x in np.ndenumerate(sfc_hght):
            # Merge the 2 meter AGL variables with the rest of the upper air profile.
            idx_aboveground = np.where(sfc_hght[index] < hght_prof)[0]
            new_hght = np.hstack((sfc_hght[index[1], index[0]]+2, hght[t, idx_aboveground, index[0], index[1]][::-1]))
            new_temp = np.hstack((sfc_temp[t, index[0], index[1]], temp[t, idx_aboveground, index[0], index[1]][::-1]))
            new_rh = np.hstack((sfc_rh[t, index[0], index[1]], rh[t, idx_aboveground, index[0], index[1]][::-1]))
            new_pres = np.hstack((sfc_pres[t, index[0], index[1]], pres[idx_aboveground][::-1]))
            new_q = utils.rh2q(new_temp, new_pres, new_rh/100.)*1000.
            
            # Store the thermodynamic profile.
            mxrs.append(new_q)
            temps.append(new_temp)
            press.append(new_pres)
            hghts.append(new_hght)
            
    # Dimensions of these Numpy arrays may be only 1-D instead of 2-D as the grid length when merging
    # the surface and upper air data together may not be equal across all grid profiles within the user-selected
    # mesh.
    temps = np.asarray(temps)
    mxrs = np.asarray(mxrs)
    press = np.asarray(press)
    hghts = np.asarray(hghts)
   
    # TODO: Figure out the format of the data that gets returned by this function.
    return temps, mxrs, press, type, path

def getARMModelPrior(model_data_path, climo_prior, yyyymmdd, size, hh, hh_delta, aeri_lon, aeri_lat):
    '''
        LEFT OVER CODE THAT CONTROLS WHICH ARM-formatted RAP/RUC FILES GET USED IN THE PRIOR GENERATON
    '''
    print "This prior is spatially centered at: " + str(aeri_lat) + ',' + str(aeri_lon)
    
    dt = datetime.strptime(yyyymmdd + str(hh), '%Y%m%d%H') - timedelta(seconds=int(hh_delta)*60*60)
    end_dt = datetime.strptime(yyyymmdd + str(hh), '%Y%m%d%H') + timedelta(seconds=int(hh_delta)*60*60)
    timed = timedelta(seconds=(60*60))
    
    print "Will be searching for model netCDF files between: " + datetime.strftime(dt, '%Y-%m-%d %H') + ' and ' + datetime.strftime(end_dt, '%Y-%m-%d %H')
    
    all_temps = None
    all_mxrs = None
    all_pres = None
    
    print "Gathering profiles within a " + str(2*size) + "x" + str(2*size) + " grid."
    types = []
    paths = []
    while dt < end_dt:
        yyyymmdd = datetime.strftime(dt, '%Y%m%d')
        hour = datetime.strftime(dt, '%H')
        print "\nGathering profiles from this date/time: " + yyyymmdd + " @ " + hour
        # THE FIRST FUNCTION ABOVE GETS CALLED HERE
        temp, mxr, pres, type, link = getARMProfiles(model_data_path, yyyymmdd, hour, aeri_lon, aeri_lat, size)
        if len(temp) == 1:
            print "\tWe weren't able to find any data from this date/time.  Let's skip it."
            dt = dt + timed
            continue
        else:
            print "\tWe were able to find data from this date/time.  Let's save it."
        paths.append(link)
        types.append(type)
        if all_temps is None:
            all_temps = temp
            all_mxrs = mxr
            #all_pres = pres
        else:
            all_temps = np.vstack((all_temps, temp))
            all_mxrs = np.vstack((all_mxrs, mxr))
            #all_pres = np.vstack((all_pres, pres))
        dt = dt + timed
    
    
    all_temps = all_temps - 273.15
    print "\nThese files were used in calculating this prior:"
    for p in paths:
        print '\t' + p
    print "\nInformation about the profiles found:"
    print "\tShape of the temperature profiles: ", all_temps.shape
    print "\tShape of the water vapor mixing ratio profiles: ", all_temps.shape
    if all_temps.shape[0] < 2000:
        print "WARNING!  THERE ARE LESS THAN 2000 PROFILES FOR THIS PRIOR."
        print "RUC/RAP ARM files are probably missing.  Canceling the generation of this prior."
        sys.exit()
    priors = np.hstack((all_temps, all_mxrs))
    
    mean = np.mean(priors, axis=0)
    print "Xa SHAPE: ", mean.shape
    cov = np.cov(priors.T)
    print "Sa SHAPE: ", cov.shape
 
    return mean, cov, climo, types, paths, yyyymmdd, hh, priors.shape[0]

# This is the code that gets called by run_prior_gen.py when we want to make a prior from
# ONLINE MOTHERLODE data. It should be used only for realtime prior generation for AERIoe.
# THIS PART OF THE CODE DOES THE REALTIME DATASET GENERATION
def getRealtimePrior(yyyymmdd, size, hh, hh_delta, aeri_lon, aeri_lat):
    '''
        Function that gets called when we want to grab data from the Motherlode site.
        used to generate semi-realtime model profiles for use in AERIoe.
    '''
    
    print "Will be pulling observations from the nearest point to: " + str(aeri_lat) + ', ' + str(aeri_lon)
    print hh_delta, yyyymmdd, hh
    dt = datetime.strptime(yyyymmdd + hh, '%Y%m%d%H') - timedelta(seconds=int(hh_delta)*60*60)

    all_temps = None
    all_mxrs = None
    all_pres = None
    print "Gathering profiles within a " + str(2*size) + "x" + str(2*size) + " grid."
    types = []
    paths = []
    dt_string = datetime.strftime(dt, '%Y%m%d%H')
    
    # Searches the MOTHERLODE UCAR THREDDS Server
    all_temps, all_mxrs, all_pres, type, link = getMotherlodeProfiles(dt_string, int(hh_delta)*2, aeri_lat, aeri_lon)
    paths = [link]
    types.append(type)

    all_temps = all_temps - 273.15
    print paths, types
    print "Shape of the temperature profiles: ", all_temps.shape
    print "Shape of the water vapor mixing ratio profiles: ", all_temps.shape
    priors = np.hstack((all_temps, all_mxrs))
    
    print "Shape of the prior: ", priors.shape
    mean = np.mean(priors, axis=0)
    print "Xa: ", mean.shape
    cov = np.cov(priors.T)
    print "Sa: ", cov.shape
    return mean, cov, climo, types, paths, yyyymmdd, hh, priors.shape[0]


def makeFile(mean, cov, dir, types, yyyymmdd, hh, paths, climo, size, t_size, n, aeri_lat, aeri_lon):
    """
        Make the netCDF Prior file!
    """
    
    priorCDF_filename = dir.strip() + '/Xa_Sa_datafile.55_levels.' + yyyymmdd + '.' + hh + '.' + types[0] + '.' + str(aeri_lat) + '.' + str(aeri_lon) + '.cdf'
    print "Saving prior file as: " + priorCDF_filename

    data = Dataset(priorCDF_filename, 'w', 'NETCDF3_CLASSIC')

    data.Date_created = datetime.strftime(datetime.now(), '%a %b %d %H:%M:%S %Y')
    data.Version = 'get_model_prior.py'
    data.Machine_used = platform.platform()
    data.model = "RUC/RAP"
    data.LBL_HOME = climo.LBL_HOME
    data.Standard_atmos = climo.Standard_atmos
    data.QC_limits_T = climo.QC_limits_T
    data.QC_limits_q = climo.QC_limits_q
    data.Comment = "Prior generated using model (" + types[0] + ") data."
    #Need to include the web links to the data used to produce these files
    #Need to include the times used to produce this file.
    data.Nsonde = str(n) + ' profiles were included in the computation of this prior dataset.'
    data.lat = aeri_lat
    data.lon = aeri_lon
    data.paths = '; '.join(paths)
    data.model_types = '; '.join(types)
    data.domain_size = str(2*size) + "x" + str(2*size)
    data.temporal_size = t_size
    data.grid_spacing = '13 km'

    print "Prior generation took place at: ", data.Date_created

    data.createDimension('height', len(mean)/2)
    data.createDimension('wnum', len(climo.variables['wnum'][:]))
    data.createDimension('height2', len(mean))

    var = data.createVariable('mean_pressure', 'f4', ('height',))
    var[:] = climo.variables['mean_pressure'][:]
    var.units = climo.variables['mean_pressure'].units
    var.long_name = climo.variables['mean_pressure'].long_name

    var = data.createVariable('height', 'f4', ('height',))
    var[:] = climo.variables['height'][:]
    var.units = climo.variables['height'].units
    var.long_name = climo.variables['height'].long_name

    var = data.createVariable('mean_temperature', 'f4', ('height',))
    var[:] = mean[:55]
    var.units = climo.variables['mean_temperature'].units
    var.long_name = climo.variables['mean_temperature'].long_name

    var = data.createVariable('mean_mixingratio', 'f4', ('height',))
    var[:] = mean[55:]
    var.units = climo.variables['mean_mixingratio'].units
    var.long_name = climo.variables['mean_mixingratio'].long_name

    var = data.createVariable('height2', 'f4', ('height2',))
    var[:] = climo.variables['height2'][:]
    var.units = climo.variables['height2'].units
    var.long_name = climo.variables['height2'].long_name

    var = data.createVariable('wnum', 'f4', ('wnum',))
    var[:] = climo.variables['wnum'][:]
    var.units = climo.variables['wnum'].units
    var.long_name = climo.variables['wnum'].long_name

    var = data.createVariable('delta_od', 'f4', ('wnum',))
    var[:] = climo.variables['delta_od'][:]
    var.units = climo.variables['delta_od'].units
    var.long_name = climo.variables['delta_od'].long_name

    var = data.createVariable('radiance_true', 'f4', ('wnum',))
    var[:] = climo.variables['radiance_true'][:]
    var.units = climo.variables['radiance_true'].units
    var.long_name = climo.variables['radiance_true'].long_name

    var = data.createVariable('radiance_fast', 'f4', ('wnum',))
    var[:] = climo.variables['radiance_fast'][:]
    var.units = climo.variables['radiance_fast'].units
    var.long_name = climo.variables['radiance_fast'].long_name

    var = data.createVariable('mean_prior', 'f4', ('height2',))
    var[:] = mean
    var.units = climo.variables['mean_prior'].units
    var.long_name = climo.variables['mean_prior'].long_name

    var = data.createVariable('covariance_prior', 'f4', ('height2','height2',))
    var[:] = cov
    var.units = climo.variables['covariance_prior'].units
    var.long_name = climo.variables['covariance_prior'].long_name

    climo.close()
    
    return priorCDF_filename

def main():
    print "END."
    
if __name__ == '__main__':
    aeri_lon = -97
    aeri_lat = 35
    print getMotherlodeProfiles('2015022412', 6, aeri_lat, aeri_lon, 10, [10,20, 30])
    #main()

