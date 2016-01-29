from netCDF4 import Dataset, date2num
#from pylab import *
import numpy as np
import sys
#from mpl_toolkits.basemap import Basemap
from datetime import datetime, timedelta
import glob
import platform
import utils

"""
    Script Name: get_model_data.py
    
    Contains two primary parsing functions:
    getARMProfiles() - which will parse out the profiles within a specified spatial domain
                       for each ARM-formatted RUC/RAP model file.
    getMotherlodeProfiles() - will parse out the profiles contained within the 13 km RAP grid
                              (used for real-time model observations for AERIoe.)
                              
    Also contains two functions that call the parsing functions and control the files that 
    get opened in order to create the model observation files for AERIoe.
    
    TODO: Include a function to get the NOMADS RAP data to generate these files for PECAN.
          Include a comment block for the other three functions in this file.
"""
 
# Height grid to interpolate the RAP/RUC profiles to:
height_grid = np.arange(0.002, 17, 0.1)

def getNOMADSRAPProfiles(yyyymmddhh, aeri_lat, aeri_lon, size):
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
    rap_path = 'http://nomads.ncdc.noaa.gov/thredds/dodsC/rap130/' + yyyymmddhh[:6] + '/' + yyyymmddhh[:8] + \
        '/rap_130_' + yyyymmddhh[:8] + '_' + yyyymmddhh[8:10] + '00_000.grb2'    
    #http://nomads.ncdc.noaa.gov/thredds/dodsC/rap130/201506/20150601/rap_130_20150601_2300_000.grb2
    # Ensure that these variables are floats and ints since they may have been passed as strings.
    aeri_lon = float(aeri_lon)
    aeri_lat = float(aeri_lat)     
    size = int(size)

    # Ensure that the file found can be opened up.
    try:
        d = Dataset(rap_path)
        # Tell the user that data has been found
        print "\tFound data: " + rap_path
    except:
        # Give messages to the user that the file couldn't be found or opened.
        # Return error codes.
        print "\tWARNING!!!"
        print "\tUnable to find data for the date and hour of: " + yyyymmddhh[:8] + ' ' + yyyymmddhh[8:10] + ' UTC'
        return None, None
        
    # Open up the 13 km latitude/longitude grid.
    ll = Dataset('13km_latlon.nc')
    lon = ll.variables['lon'][:]
    lat = ll.variables['lat'][:]
    ll.close()
    
    # Tell the user that the data is being read.
    print "\tReading in the data and converting it..."
    
    # Find the indices associated with the grid point nearest to the instrument.
    idy, idx = utils.find_index_of_nearest_xy(lon, lat, aeri_lon, aeri_lat)
   
    # Extract the pressure, temperature, RH, and height grids for the specificed spatial domain
    # as well as the surface properties of these variables.  Profiles are on an isobaric grid.
    # Variables that begin with "center_" are for the profile nearest the AERI instrument.
    #pres = d.variables['pressurepgrid'][:]
    #temp =  d.variables['tempp'][0,idy-size:idy+size,idx-size:idx+size,:]
    #center_temp =  d.variables['tempp'][0,idy,idx,:]
    #rh = d.variables['rhp'][0,idy-size:idy+size,idx-size:idx+size,:]
    #center_rh = d.variables['rhp'][0,idy,idx,:]
    #hght = d.variables['heightgpp'][0,idy-size:idy+size,idx-size:idx+size,:]
    #center_hght = d.variables['heightgpp'][0,idy,idx,:]
    
    pres = d.variables['pressure'][:]
    temp =  d.variables['Temperature'][0,:,idy-size:idy+size,idx-size:idx+size]
    center_temp =  d.variables['Temperature'][0,:,idy,idx]
    rh = d.variables['Relative_humidity'][0,:,idy-size:idy+size,idx-size:idx+size]
    center_rh = d.variables['Relative_humidity'][0,:,idy,idx]
    hght = d.variables['Geopotential_height'][0,:,idy-size:idy+size,idx-size:idx+size]
    center_hght = d.variables['Geopotential_height'][0,:,idy,idx]
    sfc_pres = d.variables['Pressure_surface'][0,idy-size:idy+size,idx-size:idx+size]
    center_sfc_pres = d.variables['Pressure_surface'][0,idy,idx]
    sfc_hght = d.variables['Geopotential_height_surface'][0,idy-size:idy+size,idx-size:idx+size]
    center_sfc_hght = d.variables['Geopotential_height_surface'][0,idy,idx]
    sfc_temp = d.variables['Temperature_height_above_ground'][0,0,idy-size:idy+size,idx-size:idx+size]
    center_sfc_temp = d.variables['Temperature_height_above_ground'][0,0,idy,idx]
    sfc_rh = d.variables['Relative_humidity_height_above_ground'][0,0,idy-size:idy+size,idx-size:idx+size]
    center_sfc_rh = d.variables['Relative_humidity_height_above_ground'][0,0,idy,idx]
    
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
        
    # Store these arrays into the distribution_profiles dictionary
    distribution_profiles = {}
    distribution_profiles['temp'] = np.asarray(temps)
    distribution_profiles['wvmr'] = np.asarray(mxrs)
    distribution_profiles['pres'] = np.asarray(press)
    distribution_profiles['hght'] = np.asarray(hghts)
    distribution_profiles['path_to_data'] = rap_path
    
    # Tell the user how many profiles were found from this:
    print "\tFound " + str(len(distribution_profiles['temp'])) + ' profiles for use.'
    
    point_profile = {}
    point_profile['temp'] = center_temp
    point_profile['wvmr'] = center_q
    point_profile['pres'] = center_pres
    point_profile['hght'] = center_hght
    point_profile['lat'] = lat[idy, idx]
    point_profile['lon'] = lon[idy, idx]
    
    return distribution_profiles, point_profile

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
        print "\t   Found data: " + files[0]
        d = Dataset(files[0])
    except:
        # Give messages to the user that the file couldn't be found or opened.
        # Return error codes.
        print "\tWARNING!!!"
        print "\tUnable to find ARM-formatted RAP/RUC dataset."
        print "\tFor the date and hour of: " + yyyymmdd + ' ' + hh + ' UTC'
        return None, None
        
    # Grab the latitude/longitude grid stored in the ARM-formatted model grid file
    lon = d.variables['longitude'][:]
    lat = d.variables['latitude'][:]
   
    # Tell the user that the data is being read.
    print "\tReading in the data and converting it..."
    
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
        
    # Store these arrays into the distribution_profiles dictionary
    distribution_profiles = {}
    distribution_profiles['temp'] = np.asarray(temps)
    distribution_profiles['wvmr'] = np.asarray(mxrs)
    distribution_profiles['pres'] = np.asarray(press)
    distribution_profiles['hght'] = np.asarray(hghts)
    distribution_profiles['path_to_data'] = files[0].split('/')[-1]
    
    # Tell the user how many profiles were found from this:
    print "\tFound " + str(len(distribution_profiles['temp'])) + ' profiles for use.'
    
    point_profile = {}
    point_profile['temp'] = center_temp
    point_profile['wvmr'] = center_q
    point_profile['pres'] = center_pres
    point_profile['hght'] = center_hght
    point_profile['lat'] = lat[idy, idx]
    point_profile['lon'] = lon[idy, idx]
    
    return distribution_profiles, point_profile

def getMotherlodeProfiles(yyyymmddhh, begin_window, end_window, aeri_lat, aeri_lon, size):
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
    try:
        d = Dataset(recent_rap_path)
        print "Found 13 km RAP data for this date on the Motherlode UCAR server."
        model_name = "RAP13km"
        path = recent_rap_path
    except:
        print "No RAP data found for this date on the Motherlode UCAR server."
        print "Data Path:", recent_rap_path
        return None, None
    
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
    temp =  d.variables['Temperature_isobaric'][begin_window:end_window,:,idy-size:idy+size,idx-size:idx+size] # 
    center_temp =  d.variables['Temperature_isobaric'][begin_window:end_window,:,idy, idx] # 
    rh = d.variables['Relative_humidity_isobaric'][begin_window:end_window,:,idy-size:idy+size,idx-size:idx+size]
    center_rh = d.variables['Relative_humidity_isobaric'][begin_window:end_window,:,idy, idx]
    hght = d.variables['Geopotential_height_isobaric'][begin_window:end_window,:,idy-size:idy+size,idx-size:idx+size]
    center_hght = d.variables['Geopotential_height_isobaric'][begin_window:end_window,:,idy ,idx]
    sfc_pres = d.variables['Pressure_surface'][begin_window:end_window, idy-size:idy+size,idx-size:idx+size]  #Pascals
    center_sfc_pres = d.variables['Pressure_surface'][begin_window:end_window, idy, idx]  #Pascals
    sfc_hght = ll.variables['Geopotential_height_surface'][idx-size:idx+size, idy-size:idy+size]
    center_sfc_hght = ll.variables['Geopotential_height_surface'][idx, idy]
    #sfc_hght = d.variables['Geopotential_height_surface'][:time_window, idy-size:idy+size,idx-size:idx+size]
    sfc_temp = d.variables['Temperature_height_above_ground'][begin_window:end_window,0, idy-size:idy+size,idx-size:idx+size]
    center_sfc_temp = d.variables['Temperature_height_above_ground'][begin_window:end_window,0, idy, idx]
    sfc_rh = d.variables['Relative_humidity_height_above_ground'][begin_window:end_window,0, idy-size:idy+size,idx-size:idx+size]
    center_sfc_rh = d.variables['Relative_humidity_height_above_ground'][begin_window:end_window,0, idy, idx]
   
    # Merge the model profile data together to get the thermodynamic profile nearest to the AERI location.
    hght_prof = center_hght[0, :]
    idx_aboveground = np.where(pres < center_sfc_pres[0])[0]
    center_hght = (np.hstack((center_sfc_hght+2, center_hght[0, idx_aboveground][::-1])) - center_sfc_hght+2)/1000.
    center_temp = np.hstack((center_sfc_temp[0], center_temp[0, idx_aboveground][::-1]))
    center_rh = np.hstack((center_sfc_rh[0], center_rh[0, idx_aboveground][::-1]))
    center_pres = np.hstack((center_sfc_pres[0], pres[idx_aboveground][::-1]))
    center_q = utils.rh2q(center_temp, center_pres, center_rh/100.)*1000.
    
    # Close the 13 km RAP dataset.
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
            new_hght = (np.hstack((sfc_hght[index[1], index[0]]+2, hght[t, idx_aboveground, index[0], index[1]][::-1])) - sfc_hght[index[1], index[0]]+2)/1000.
            new_temp = np.hstack((sfc_temp[t, index[0], index[1]], temp[t, idx_aboveground, index[0], index[1]][::-1]))
            new_rh = np.hstack((sfc_rh[t, index[0], index[1]], rh[t, idx_aboveground, index[0], index[1]][::-1]))
            new_pres = np.hstack((sfc_pres[t, index[0], index[1]], pres[idx_aboveground][::-1]))
            new_q = utils.rh2q(new_temp, new_pres, new_rh/100.)*1000.
            
            # Store the thermodynamic profile.
            mxrs.append(new_q)
            temps.append(new_temp)
            press.append(new_pres)
            hghts.append(new_hght)
            
    # Store these arrays into the distribution_profiles dictionary
    distribution_profiles = {}
    distribution_profiles['temp'] = np.asarray(temps)
    distribution_profiles['wvmr'] = np.asarray(mxrs)
    distribution_profiles['pres'] = np.asarray(press)
    distribution_profiles['hght'] = np.asarray(hghts)
    distribution_profiles['path_to_data'] = recent_rap_path
    
    point_profile = {}
    point_profile['temp'] = center_temp
    point_profile['wvmr'] = center_q
    point_profile['pres'] = center_pres
    point_profile['hght'] = center_hght
    point_profile['lat'] = ll.variables['lat'][idy, idx]
    point_profile['lon'] = ll.variables['lon'][idy, idx]
    
    # Close the terrain dataset.
    ll.close()
    
    return distribution_profiles, point_profile

def getNOMADSModelObs(begin_dt, end_dt, temporal_mesh_size, spatial_mesh_size, aeri_lon, aeri_lat):
    '''
        LEFT OVER CODE THAT CONTROLS WHICH ARM-formatted RAP/RUC FILES GET USED IN THE PRIOR GENERATON
    '''
    print "This model sounding is spatially centered at: " + str(aeri_lat) + ',' + str(aeri_lon)
    delta = timedelta(seconds=60*60) # Hour delta used to iterate throughout the files
    
    # Tell the user what the range of data the program will look for.
    print "Will be searching for model netCDF files between: " + datetime.strftime(begin_dt, '%Y-%m-%d %H') + ' and ' + datetime.strftime(end_dt, '%Y-%m-%d %H')
    print "Gathering profiles within a " + str(2*spatial_mesh_size) + "x" + str(2*spatial_mesh_size) + " grid."

    # Determine the temporal bounds of the data we need to collect.
    lower_bound_dt = begin_dt - timedelta(seconds=60*60*temporal_mesh_size)
    upper_bound_dt = end_dt + timedelta(seconds=60*60*temporal_mesh_size)
    
    # Arrays to store dictionaries 
    num_times = (upper_bound_dt - lower_bound_dt).seconds * (1./3600.) + (upper_bound_dt - lower_bound_dt).days * 24
    dists = np.empty(num_times+1, dtype=dict)
    points = np.empty(num_times+1, dtype=dict)
    dts = np.empty(num_times+1, dtype=object)
    
    # Begin looping over the time frame the observation files encompass.
    # Save the data into the numpy arrays.
    cur_dt = lower_bound_dt
    count = 0
    
    paths = []
    while cur_dt <= upper_bound_dt :
        print "\nGathering profiles from this date/time: " + datetime.strftime(cur_dt, '%Y%m%d %H UTC')
        #dist, point = getARMProfiles(model_data_path, , datetime.strftime(cur_dt,'%H'),aeri_lon, aeri_lat, )
        dist, point = getNOMADSRAPProfiles(datetime.strftime(cur_dt, '%Y%m%d%H'), aeri_lat, aeri_lon, spatial_mesh_size)
        paths.append(dist['path_to_data'])
        dists[count] = dist
        points[count] = point
        dts[count] = cur_dt
        cur_dt = cur_dt + delta
        count = count + 1
    
    paths = ', '.join(paths)
    
    temperature = np.zeros((len(dts[temporal_mesh_size:len(dts)-temporal_mesh_size]), len(height_grid)))
    wvmr = np.zeros((len(dts[temporal_mesh_size:len(dts)-temporal_mesh_size]), len(height_grid)))
    pressure = np.zeros((len(dts[temporal_mesh_size:len(dts)-temporal_mesh_size]), len(height_grid)))
    temperature_sigma = np.zeros((len(dts[temporal_mesh_size:len(dts)-temporal_mesh_size]), len(height_grid)))
    wvmr_sigma = np.zeros((len(dts[temporal_mesh_size:len(dts)-temporal_mesh_size]), len(height_grid)))
    
    output = {}
    # Loop over the timeframe.
    for i in np.arange(temporal_mesh_size, len(dts) - temporal_mesh_size, 1):
        index_range = np.arange(i - temporal_mesh_size, i+temporal_mesh_size+1, 1)
        
        # Try to save the interpolated temperature, water vapor mixing ratio, and pressure profiles
        try:
            temperature[i-temporal_mesh_size,:] = np.interp(height_grid, points[i]['hght'], points[i]['temp'])
            wvmr[i-temporal_mesh_size,:] = np.interp(height_grid, points[i]['hght'], points[i]['wvmr'])
            pressure[i-temporal_mesh_size,:] = np.interp(height_grid, points[i]['hght'], points[i]['pres'])
        except Exception,e:
            # If there's an issue with loading in the data for this time, then this exception will catch it.
            # this will skip calculating the standard deviation too, because we won't need that.
            print e
            continue
        
        # Pull out the latitude and longitude point.
        lat = points[i]['lat']
        lon = points[i]['lon']
        
        temp_dist = []
        wvmr_dist = []
        # Loop over the temporal window we'll be using to calculate the standard deviation.
        for j in index_range:
            # Loop over the spatial distribution of profiles for time index j.
            for k in range(len(dists[i]['temp'])):
                # Interpolate the profiles and save to the array used in calculating the standard deviation.
                temp_dist.append(np.interp(height_grid, dists[i]['hght'][k], dists[i]['temp'][k]))
                wvmr_dist.append(np.interp(height_grid, dists[i]['hght'][k], dists[i]['wvmr'][k]))
        # Calculate the standard deviation profiles for temperature and water vapor mixing ratio and save it.
        temp_std = np.std(np.asarray(temp_dist), axis=0)
        wvmr_std = np.std(np.asarray(wvmr_dist), axis=0)
        temperature_sigma[i-temporal_mesh_size,:] = temp_std
        wvmr_sigma[i-temporal_mesh_size,:] = wvmr_std
        
    output['pressure'] = np.ma.masked_where(pressure == 0, pressure)
    output['temperature'] = np.ma.masked_where(np.ma.getmask(output['pressure']), temperature)
    output['wvmr'] = np.ma.masked_where(np.ma.getmask(output['pressure']), wvmr)
    output['temperature_sigma'] = np.ma.masked_where(np.ma.getmask(output['pressure']), temperature_sigma)
    output['wvmr_sigma'] = np.ma.masked_where(np.ma.getmask(output['pressure']), wvmr_sigma)
    output['height'] = height_grid
    output['paths_to_data'] = paths
    output['gridpoint_lat'] = lat
    output['gridpoint_lon'] = lon
    output['data_type'] = np.ones(len(temperature))
    output['dts'] = dts[temporal_mesh_size:len(dts)-temporal_mesh_size]
    
    return output
    
    
def getARMModelObs(model_data_path, begin_dt, end_dt, temporal_mesh_size, spatial_mesh_size, aeri_lon, aeri_lat):
    '''
        LEFT OVER CODE THAT CONTROLS WHICH ARM-formatted RAP/RUC FILES GET USED IN THE PRIOR GENERATON
    '''
    print "This model sounding is spatially centered at: " + str(aeri_lat) + ',' + str(aeri_lon)
    delta = timedelta(seconds=60*60) # Hour delta used to iterate throughout the files
    
    # Tell the user what the range of data the program will look for.
    print "Will be searching for model netCDF files between: " + datetime.strftime(begin_dt, '%Y-%m-%d %H') + ' and ' + datetime.strftime(end_dt, '%Y-%m-%d %H')
    print "Gathering profiles within a " + str(2*spatial_mesh_size) + "x" + str(2*spatial_mesh_size) + " grid."

    # Determine the temporal bounds of the data we need to collect.
    lower_bound_dt = begin_dt - timedelta(seconds=60*60*temporal_mesh_size)
    upper_bound_dt = end_dt + timedelta(seconds=60*60*temporal_mesh_size)
    
    # Arrays to store dictionaries 
    num_times = (upper_bound_dt - lower_bound_dt).seconds * (1./3600.) + (upper_bound_dt - lower_bound_dt).days * 24
    dists = np.empty(num_times+1, dtype=dict)
    points = np.empty(num_times+1, dtype=dict)
    dts = np.empty(num_times+1, dtype=object)
    
    # Begin looping over the time frame the observation files encompass.
    # Save the data into the numpy arrays.
    cur_dt = lower_bound_dt
    count = 0
    
    paths = []
    while cur_dt <= upper_bound_dt :
        print "\nGathering profiles from this date/time: " + datetime.strftime(cur_dt, '%Y%m%d %H UTC')
        dist, point = getARMProfiles(model_data_path, datetime.strftime(cur_dt, '%Y%m%d'), datetime.strftime(cur_dt,'%H'),aeri_lon, aeri_lat, spatial_mesh_size)
        if dist is not None:
            paths.append(dist['path_to_data'])
        dists[count] = dist
        points[count] = point
        dts[count] = cur_dt
        cur_dt = cur_dt + delta
        count = count + 1
    
    paths = ', '.join(paths)
    
    temperature = np.zeros((len(dts[temporal_mesh_size:len(dts)-temporal_mesh_size]), len(height_grid)))
    wvmr = np.zeros((len(dts[temporal_mesh_size:len(dts)-temporal_mesh_size]), len(height_grid)))
    pressure = np.zeros((len(dts[temporal_mesh_size:len(dts)-temporal_mesh_size]), len(height_grid)))
    temperature_sigma = np.zeros((len(dts[temporal_mesh_size:len(dts)-temporal_mesh_size]), len(height_grid)))
    wvmr_sigma = np.zeros((len(dts[temporal_mesh_size:len(dts)-temporal_mesh_size]), len(height_grid)))
    
    output = {}
    # Loop over the timeframe.
    for i in np.arange(temporal_mesh_size, len(dts) - temporal_mesh_size, 1):
        index_range = np.arange(i - temporal_mesh_size, i+temporal_mesh_size+1, 1)
        
        # Try to save the interpolated temperature, water vapor mixing ratio, and pressure profiles
        try:
            temperature[i-temporal_mesh_size,:] = np.interp(height_grid, points[i]['hght'], points[i]['temp'])
            wvmr[i-temporal_mesh_size,:] = np.interp(height_grid, points[i]['hght'], points[i]['wvmr'])
            pressure[i-temporal_mesh_size,:] = np.interp(height_grid, points[i]['hght'], points[i]['pres'])
        except Exception,e:
            # If there's an issue with loading in the data for this time, then this exception will catch it.
            # this will skip calculating the standard deviation too, because we won't need that.
            continue
        
        # Pull out the latitude and longitude point.
        lat = points[i]['lat']
        lon = points[i]['lon']
        
        temp_dist = []
        wvmr_dist = []
        # Loop over the temporal window we'll be using to calculate the standard deviation.
        for j in index_range:
            # Loop over the spatial distribution of profiles for time index j.
            for k in range(len(dists[i]['temp'])):
                # Interpolate the profiles and save to the array used in calculating the standard deviation.
                temp_dist.append(np.interp(height_grid, dists[i]['hght'][k], dists[i]['temp'][k]))
                wvmr_dist.append(np.interp(height_grid, dists[i]['hght'][k], dists[i]['wvmr'][k]))
        # Calculate the standard deviation profiles for temperature and water vapor mixing ratio and save it.
        temp_std = np.std(np.asarray(temp_dist), axis=0)
        wvmr_std = np.std(np.asarray(wvmr_dist), axis=0)
        temperature_sigma[i-temporal_mesh_size,:] = temp_std
        wvmr_sigma[i-temporal_mesh_size,:] = wvmr_std
        
    output['pressure'] = np.ma.masked_where(pressure == 0, pressure)
    output['temperature'] = np.ma.masked_where(np.ma.getmask(output['pressure']), temperature)
    output['wvmr'] = np.ma.masked_where(np.ma.getmask(output['pressure']), wvmr)
    output['temperature_sigma'] = np.ma.masked_where(np.ma.getmask(output['pressure']), temperature_sigma)
    output['wvmr_sigma'] = np.ma.masked_where(np.ma.getmask(output['pressure']), wvmr_sigma)
    output['height'] = height_grid
    output['paths_to_data'] = paths
    output['gridpoint_lat'] = lat
    output['gridpoint_lon'] = lon
    output['data_type'] = np.ones(len(temperature))
    output['dts'] = dts[temporal_mesh_size:len(dts)-temporal_mesh_size]
    
    total_missing = np.ma.count_masked(output['pressure'], axis=0)
    total = np.ma.count(output['pressure'], axis=0)
    
    if total_missing[0] != 0:
        print "WARNING!  Some files were missing when converting the data.  A total of " + str(total_missing[0]) + ' profiles out of ' + str(len(output['dts'])) + ' are missing!'
    
    return output

# This is the code that gets called by run_prior_gen.py when we want to make a prior from
# ONLINE MOTHERLODE data. It should be used only for realtime prior generation for AERIoe.
# THIS PART OF THE CODE DOES THE REALTIME DATASET GENERATION
def getRealtimeProfiles(begin_dt, end_dt, temporal_mesh_size, spatial_mesh_size, aeri_lon, aeri_lat):
    '''
        Function that gets called when we want to grab data from the Motherlode site.
        used to generate semi-realtime model profiles for use in AERIoe.
    '''
    print "This model sounding is spatially centered at: " + str(aeri_lat) + ',' + str(aeri_lon)
    delta = timedelta(seconds=60*60) # Hour delta used to iterate throughout the files
    
    # Tell the user what the range of data the program will look for.
    print "Will be searching for model netCDF files between: " + datetime.strftime(begin_dt, '%Y-%m-%d %H') + ' and ' + datetime.strftime(end_dt, '%Y-%m-%d %H')
    print "Gathering profiles within a " + str(2*spatial_mesh_size) + "x" + str(2*spatial_mesh_size) + " grid."

    # Determine the temporal bounds of the data we need to collect.
    lower_bound_dt = begin_dt - timedelta(seconds=60*60*temporal_mesh_size)
    upper_bound_dt = end_dt + timedelta(seconds=60*60*temporal_mesh_size)
    
    # Arrays to store dictionaries 
    num_times = (upper_bound_dt - lower_bound_dt).seconds * (1./3600.) + (upper_bound_dt - lower_bound_dt).days * 24
    dists = np.empty(num_times+24, dtype=dict)
    points = np.empty(num_times+24, dtype=dict)
    dts = []
    
    # Begin looping over the time frame the observation files encompass.
    # Save the data into the numpy arrays.
    cur_dt = lower_bound_dt
    count = 0
    
    paths = []
    data_type = []
    use_forecast = False
    # Get all of the analysis profiles
    while cur_dt <= upper_bound_dt :
        print "\nGathering profiles from this date/time: " + datetime.strftime(cur_dt, '%Y%m%d %H UTC')
        dist, point = getMotherlodeProfiles(datetime.strftime(cur_dt, '%Y%m%d%H'), 0, 1, aeri_lat, aeri_lon, spatial_mesh_size)
        if dist == None:
            # This means that a file couldn't be found for this time.
            # This means that we should break the loop because no data for this hour is available.
            # This also means that we should use forecast data to fill in the distribution
            print "Unable to find:", cur_dt
            use_forecast = True
            break
        data_type.append(1) # Means that this an analysis
        paths.append(dist['path_to_data'])
        dists[count] = dist
        points[count] = point
        dts.append(cur_dt)
        cur_dt = cur_dt + delta
        count = count + 1
    
    if use_forecast == True:
        # Go back to the last file that actually existed and had data
        cur_dt = cur_dt - delta 

        # Get all of the forecast profiles from that file.
        for i in range(temporal_mesh_size):
            dist, point = getMotherlodeProfiles(datetime.strftime(cur_dt, '%Y%m%d%H'), i, i+1, aeri_lat, aeri_lon, spatial_mesh_size)
            data_type.append(2) # Means that this a forecast value
            paths.append(dist['path_to_data'])
            dists[count+i] = dist
            points[count+i] = point
            dts.append(cur_dt+delta)
            count = count + 1
            
    # Filter out the array elements that were never filled with profile information.
    idx_filter = [i for i, item in enumerate(points) if item is not None]
    points = points[idx_filter]
    dists = dists[idx_filter]

    paths = ', '.join(paths)
    
    # Set up final data arrays to be saved to the file.
    temperature = np.zeros((len(dts[temporal_mesh_size:len(dts)-temporal_mesh_size]), len(height_grid)))
    wvmr = np.zeros((len(dts[temporal_mesh_size:len(dts)-temporal_mesh_size]), len(height_grid)))
    pressure = np.zeros((len(dts[temporal_mesh_size:len(dts)-temporal_mesh_size]), len(height_grid)))
    temperature_sigma = np.zeros((len(dts[temporal_mesh_size:len(dts)-temporal_mesh_size]), len(height_grid)))
    wvmr_sigma = np.zeros((len(dts[temporal_mesh_size:len(dts)-temporal_mesh_size]), len(height_grid)))
    
    output = {}
    # Loop over the timeframe.
    for i in np.arange(temporal_mesh_size, len(dts) - temporal_mesh_size, 1):
        index_range = np.arange(i - temporal_mesh_size, i+temporal_mesh_size+1, 1)
        
        # Try to save the interpolated temperature, water vapor mixing ratio, and pressure profiles
        try:
            temperature[i-temporal_mesh_size,:] = np.interp(height_grid, points[i]['hght'], points[i]['temp'])
            wvmr[i-temporal_mesh_size,:] = np.interp(height_grid, points[i]['hght'], points[i]['wvmr'])
            pressure[i-temporal_mesh_size,:] = np.interp(height_grid, points[i]['hght'], points[i]['pres'])
        except Exception,e:
            # If there's an issue with loading in the data for this time, then this exception will catch it.
            # this will skip calculating the standard deviation too, because we won't need that.
            print e
            continue
        
        # Pull out the latitude and longitude point.
        lat = points[i]['lat']
        lon = points[i]['lon']
        
        temp_dist = []
        wvmr_dist = []
        # Loop over the temporal window we'll be using to calculate the standard deviation.
        for j in index_range:
            # Loop over the spatial distribution of profiles for time index j.
            for k in range(len(dists[i]['temp'])):
                # Interpolate the profiles and save to the array used in calculating the standard deviation.
                temp_dist.append(np.interp(height_grid, dists[i]['hght'][k], dists[i]['temp'][k]))
                wvmr_dist.append(np.interp(height_grid, dists[i]['hght'][k], dists[i]['wvmr'][k]))
        # Calculate the standard deviation profiles for temperature and water vapor mixing ratio and save it.
        temp_std = np.std(np.asarray(temp_dist), axis=0)
        wvmr_std = np.std(np.asarray(wvmr_dist), axis=0)
        temperature_sigma[i-temporal_mesh_size,:] = temp_std
        wvmr_sigma[i-temporal_mesh_size,:] = wvmr_std
    
    output['pressure'] = np.ma.masked_where(pressure == 0, pressure)
    output['temperature'] = np.ma.masked_where(np.ma.getmask(output['pressure']), temperature)
    output['wvmr'] = np.ma.masked_where(np.ma.getmask(output['pressure']), wvmr)
    output['temperature_sigma'] = np.ma.masked_where(np.ma.getmask(output['pressure']), temperature_sigma)
    output['wvmr_sigma'] = np.ma.masked_where(np.ma.getmask(output['pressure']), wvmr_sigma)
    output['height'] = height_grid
    output['paths_to_data'] = paths
    output['gridpoint_lat'] = lat
    output['gridpoint_lon'] = lon
    output['data_type'] = np.ones(len(temperature))
    output['dts'] = dts[temporal_mesh_size:len(dts)-temporal_mesh_size]
   
    return output


def makeFile(output):
    """
        Make the netCDF file containing the model observations!
    """
    epoch_time = date2num(output['dts'], 'seconds since 1979-01-01 00:00:00+00:00')
    priorCDF_filename = output['output_dir'] + '/RRmodelsoundings.' + datetime.strftime(output['dts'][0], '%Y%m%d') + '.' + datetime.strftime(output['dts'][0], '%H') + '.' + str(output['aeri_lat']) + '.' + str(output['aeri_lon']) + '.cdf'
    print "Saving prior file as: " + priorCDF_filename

    data = Dataset(priorCDF_filename, 'w', 'NETCDF3_CLASSIC')

    data.Date_created = datetime.strftime(datetime.now(), '%a %b %d %H:%M:%S %Y')
    data.Machine_used = platform.platform()
    data.model = "RUC/RAP"
    #Need to include the web links to the data used to produce these files
    #Need to include the times used to produce this file.
    data.aeri_lat = output['aeri_lat']
    data.aeri_lon = output['aeri_lon']
    data.gridpoint_lat = output['gridpoint_lat']
    data.gridpoint_lon = output['gridpoint_lon']
    data.paths = output['paths_to_data']
    data.arm_model_dir = output['arm_model_dir']
    data.domain_size = str(2*output['spatial_mesh_size']) + "x" + str(2*output['spatial_mesh_size'])
    data.temporal_size = output['temporal_mesh_size']

    print "Prior generation took place at: ", data.Date_created

    data.createDimension('time', len(output['temperature']))
    data.createDimension('height', len(output['height']))

    var = data.createVariable('base_time', 'i4')
    var[:] = epoch_time[0]
    var.units = 'seconds since 1970-01-01 00:00:00+00:00'
    var.long_name = 'epoch time'

    var = data.createVariable('time_offset', 'f4', ('time',))
    var[:] = epoch_time - epoch_time[0]
    var.units = 's'
    var.long_name = 'Time offset from base_time'
    
    var = data.createVariable('height', 'f4', ('height',))
    var[:] = output['height']
    var.units = 'km AGL'
    var.long_name = 'height grid of the sounding'

    var = data.createVariable('temperature', 'f4', ('time','height',))
    var[:] = output['temperature'] - 273.15
    var.units = 'C'
    var.long_name = "temperature"

    var = data.createVariable('waterVapor', 'f4', ('time','height',))
    var[:] = output['wvmr']
    var.units = 'g/kg'
    var.long_name = 'water vapor mixing ratio'

    var = data.createVariable('pressure', 'f4', ('time','height',))
    var[:] = output['pressure']
    var.units = 'mb'
    var.long_name = 'air pressure'

    var = data.createVariable('sigma_temperature', 'f4', ('time','height',))
    var[:] = output['temperature_sigma']
    var.units = 'C'
    var.long_name = '1-sigma uncertainty in temperature'

    var = data.createVariable('sigma_waterVapor', 'f4', ('time','height',))
    var[:] = output['wvmr_sigma']
    var.units = 'g/kg'
    var.long_name = '1-sigma uncertainty in water vapor mixing ratio'

    var = data.createVariable('lat', 'f4')
    var[:] = output['gridpoint_lat']
    var.units = 'degrees north'
    var.long_name = 'latitude'

    var = data.createVariable('lon', 'f4')
    var[:] = output['gridpoint_lon']
    var.units = 'degrees east'
    var.long_name = 'longitude'
    
    data.close()
    
    return priorCDF_filename

def main():
    print "END."
    
if __name__ == '__main__':
    aeri_lon = -97
    aeri_lat = 35
    print getMotherlodeProfiles('2015022412', 6, aeri_lat, aeri_lon, 10, [10,20, 30])
    #main()

