from netCDF4 import Dataset
from pylab import *
import numpy as np
import sys
#from mpl_toolkits.basemap import Basemap
from datetime import datetime, timedelta
import glob
import platform

"""
    Script Name: get_model_prior.py
    
"""

def getARMProfiles(rap_path, yyyymmdd, hh,  aeri_lon, aeri_lat, size, aerioe_hght):
    # THIS FUNCTION GETS CALLED BY getARMModelPrior() [this is defined below]
    # This is the function that is called to pull out the profiles from the 
    # ARM RUC/RAP "syn" files contained in the rap_path directory.
    # This code was moved here from the file get_ruca_prior.py on 3/26/2014
    #print rap_path + "/*syn*" + yyyymmdd + '.' + hh + '*.cdf'
    print "\tSearching for \"syn\" files."
    files = glob.glob(rap_path.strip() + '/*syn*' + yyyymmdd + '.' + hh + '*.cdf') 
    if len(files) == 0:
        print "\tNone, searching for \"all\" files."
        files = glob.glob(rap_path.strip() + '/*all*' + yyyymmdd + '.' + hh + '*.cdf') 
    
    aeri_lon = float(aeri_lon)
    aeri_lat = float(aeri_lat)     
    size = int(size)

    try:
        d = Dataset(files[0])
    except:
        print "\t\tUnable to find ARM data to generate prior."
        print "\t\tARM RAP/RUC \"syn\" data needs to be in the directory: " + rap_path
        print "\t\tFor the date and hour of: " + yyyymmdd + ' ' + hh + ' UTC'
        return [-9999],[-9999], [-9999], [-9999], [-9999]
        #sys.exit()
    print "\t\tFOUND DATA FOR PRIOR: " + files[0]
        
    lon = d.variables['longitude'][:]
    lat = d.variables['latitude'][:]
    #print d.variables.keys()
    idy, idx = find_index_of_nearest_xy(lon, lat, aeri_lon, aeri_lat)
    
    pres = d.variables['pressurepgrid'][:]
    temp =  d.variables['tempp'][0,idy-size:idy+size,idx-size:idx+size,:]
    rh = d.variables['rhp'][0,idy-size:idy+size,idx-size:idx+size,:]
    hght = d.variables['heightgpp'][0,idy-size:idy+size,idx-size:idx+size,:]
    sfc_pres = d.variables['pressuresrf'][0,idy-size:idy+size,idx-size:idx+size]/100.
    sfc_hght = d.variables['heightsrf'][idy-size:idy+size,idx-size:idx+size]
    sfc_temp = d.variables['temp2m'][0,idy-size:idy+size,idx-size:idx+size]
    sfc_rh = d.variables['rh2m'][0,idy-size:idy+size,idx-size:idx+size]
    
    d.close()
    mxrs = []
    temps = []
    press = []
    #print "sfc_hght.shape: ", sfc_hght.shape
    for index, x in np.ndenumerate(sfc_hght):
        idx_aboveground = np.where(sfc_hght[index] < hght[index[0], index[1],:])[0]
 
        new_hght = (np.hstack((sfc_hght[index]+2, hght[index[0], index[1], idx_aboveground])) - sfc_hght[index])/1000.
        new_temp = np.hstack((sfc_temp[index], temp[index[0], index[1], idx_aboveground]))
        new_rh = np.hstack((sfc_rh[index], rh[index[0], index[1], idx_aboveground]))
        new_pres = np.hstack((sfc_pres[index], pres[idx_aboveground]))
        new_q = rh2q(new_temp, new_pres*100., new_rh/100.)*1000.
        oe_mxr = np.interp(aerioe_hght, new_hght, new_q)
        oe_temp = np.interp(aerioe_hght, new_hght, new_temp)
        oe_pres = np.interp(aerioe_hght, new_hght, new_pres)
        mxrs.append(oe_mxr)
        temps.append(oe_temp)
        press.append(new_pres)
    temps = np.asarray(temps)
    mxrs = np.asarray(mxrs)
    type = files[0].split('/')[-1].split('.')[0]   
    #print temps.shape, mxrs.shape, np.asarray(press).shape, type, files[0]
    
    return temps, mxrs, press, type, files[0]

def getMotherlodeProfiles(yyyymmddhh, time_window, aeri_lat, aeri_lon):
    #This function is used to check for recent RAP data on the Motherlode server
    #This is the function that should be called if we need to run AERIoe in real-time.
    #
    # THIS FUNCTION IS WORKING ON 2/24/2015
    # GETS WHEN I WANT REALTIME RAP OBSERVATIONS
    # 
    # FORECASTS AND ANALYSES HAVE A 2 WEEK LIFETIME ON THIS SERVER

    recent_rap_path = 'http://thredds.ucar.edu/thredds/dodsC/grib/NCEP/RAP/CONUS_13km/RR_CONUS_13km_' + \
        yyyymmddhh[:8] + '_' + yyyymmddhh[8:10] + '00.grib2/GC'
    print recent_rap_path
    try:
        d = Dataset(recent_rap_path)
        print d.variables.keys()
        print "Found 13 km RAP data for this date on the Motherlode UCAR server."
        type = "RAP13km"
        path = recent_rap_path
    except:
        print "No RAP data found for this date on the Motherlode UCAR server."
        return None
    
    ll = Dataset('13km_latlon.nc')
    lon = ll.variables['lon'][:]
    lat = ll.variables['lat'][:]
    aeri_lon = float(aeri_lon)
    aeri_lat = float(aeri_lat)
    idy, idx = find_index_of_nearest_xy(lon, lat, aeri_lon, aeri_lat)

    pres = d.variables['isobaric'][:] #Pascals
    temp =  d.variables['Temperature_isobaric'][:time_window,:,idy-size:idy+size,idx-size:idx+size] # 
    center_temp =  d.variables['Temperature_isobaric'][:time_window,:,idy, idx] # 
    rh = d.variables['Relative_humidity_isobaric'][:time_window,:,idy-size:idy+size,idx-size:idx+size]
    center_rh = d.variables['Relative_humidity_isobaric'][:time_window,:,idy, idx]
    hght = d.variables['Geopotential_height_isobaric'][:time_window,:,idy-size:idy+size,idx-size:idx+size]
    center_hght = d.variables['Geopotential_height_isobaric'][:time_window,:,idy ,idx]
    sfc_pres = d.variables['Pressure_surface'][:time_window, idy-size:idy+size,idx-size:idx+size]  #Pascals
    center_sfc_pres = d.variables['Pressure_surface'][:time_window, idy, idx]  #Pascals
    sfc_hght = ll.variables['Geopotential_height_surface'][idy-size:idy+size, idx-size:idx+size]
    center_sfc_hght = ll.variables['Geopotential_height_surface'][idx, idy]
    #sfc_hght = d.variables['Geopotential_height_surface'][:time_window, idy-size:idy+size,idx-size:idx+size]
    sfc_temp = d.variables['Temperature_height_above_ground'][:time_window,0, idy-size:idy+size,idx-size:idx+size]
    center_sfc_temp = d.variables['Temperature_height_above_ground'][:time_window,0, idy, idx]
    sfc_rh = d.variables['Relative_humidity_height_above_ground'][:time_window,0, idy-size:idy+size,idx-size:idx+size]
    center_sfc_rh = d.variables['Relative_humidity_height_above_ground'][:time_window,0, idy, idx]
    
    print "i, TEMP, RH, HGHT, PRES"
    for i in range(len(center_temp[0])):
        print i, center_temp[0][i], center_rh[0][i], center_hght[0][i], pres[i]
    print center_sfc_pres[0], center_sfc_temp[0], center_sfc_hght, center_sfc_rh[0]
    
    hght_prof = center_hght[0, :]
    idx_aboveground = np.where(pres < center_sfc_pres[0])[0]
    print idx_aboveground
    new_hght = np.hstack((center_sfc_hght+2, center_hght[0, idx_aboveground][::-1]))
    new_temp = np.hstack((center_sfc_temp[0], center_temp[0, idx_aboveground][::-1]))
    new_rh = np.hstack((center_sfc_rh[0], center_rh[0, idx_aboveground][::-1]))
    new_pres = np.hstack((center_sfc_pres[0], pres[idx_aboveground][::-1]))
    new_q = rh2q(new_temp, new_pres, new_rh/100.)*1000.
    
    print "i, TEMP, RH, HGHT, PRES, WVMR"
    for i in range(len(new_q)):
        print i, new_temp[i], new_rh[i], new_hght[i], new_pres[i], new_q[i]

    stop
    ll.close()
    d.close()
    
    mxrs = []
    temps = []
    press = []
    for t in range(len(temp)):
        for index, x in np.ndenumerate(sfc_hght):
            hght_prof = hght[t, :, index[0], index[1]]
            idx_aboveground = np.where(sfc_hght[index] < hght_prof)[0]
            new_hght = (np.hstack((sfc_hght[index]+2, hght[t, idx_aboveground, index[0], index[1]][::-1])) - sfc_hght[index])/1000.
            new_temp = np.hstack((sfc_temp[t, index[0], index[1]], temp[t, idx_aboveground, index[0], index[1]][::-1]))
            new_rh = np.hstack((sfc_rh[t, index[0], index[1]], rh[t, idx_aboveground, index[0], index[1]][::-1]))
            new_pres = np.hstack((sfc_pres[t, index[0], index[1]], pres[idx_aboveground][::-1]))
            new_q = rh2q(new_temp, new_pres, new_rh/100.)*1000.
            print new_q, new_temp, new_hght, new_pres 
            print new_q.shape, new_temp.shape, new_hght.shape, new_pres.shape

            stop 
            oe_mxr = np.interp(aerioe_hght, new_hght, new_q)
            oe_temp = np.interp(aerioe_hght, new_hght, new_temp)
            oe_pres = np.interp(aerioe_hght, new_hght, new_pres)
            mxrs.append(oe_mxr)
            temps.append(oe_temp)
            press.append(new_pres)
    
    temps = np.asarray(temps)
    mxrs = np.asarray(mxrs)
    
    return temps, mxrs, press, type, path
"""
#This function searches for the RAP/RUC forecast data on the NOMADS THREDDS SERVER and builds a prior from it.
def getHourlyProfiles(yyyymmddhh, fff, aeri_lat, aeri_lon, size, aerioe_hght):
    #This function is used if the RAP data cannot be found on the Motherlode server.

    rap_path = 'http://nomads.ncdc.noaa.gov/thredds/dodsC/rap130/' + yyyymmddhh[:6] + '/' + yyyymmddhh[:8] + \
        '/rap_130_' + yyyymmddhh[:8] + '_' + yyyymmddhh[8:10] + '00_' + fff + '.grb2'
   
    #ruc_path = 'http://nomads.ncdc.noaa.gov/thredds/dodsC/ruc252/' + yyyymmddhh[:6] + '/' + yyyymmddhh[:8] + \
    #    '/ruc2_252_' + yyyymmddhh[:8] + '_' + yyyymmddhh[8:10] + '00_' + fff + '.grb'
    
    found = False
    for i in range(1,6):
        try:
            print "Attempting to find: " + rap_path
            d = Dataset(rap_path)
            print "RAP 13 km data found."
            type = "RAP13km"
            path = rap_path
            latlon_path = '13km_latlon.nc'
            found = True
            break
        except:
            print "Unable to find data on attempt " + str(i)
            continue
    if found == False:
        for i in range(1,6):
            try:
                print "RAP data not found."
                print "Attempting to find: " + ruc_path
                d = Dataset(ruc_path)
                print "RUC 20 km data found."
                type = "RUC20km"
                path = ruc_path
                latlon_path = '20km_latlon.nc'
                found = True
                break
            except:
                print "Unable to find data on attempt " + str(i)
                continue

    if found == False:
        print "Data not found on the server, aborting."
        print yyyymmddhh + 'F' + fff
        print "\n"
        print "Maybe you should try using a GFS/CFS based prior instead?"
        sys.exit()

    try:
        pres = d.variables['pressure'][:]
    except:
        pres = d.variables['isobaric'][:]*100.
    
    ll = Dataset(latlon_path)
    lon = ll.variables['lon'][:]
    lat = ll.variables['lat'][:]
    idy, idx = find_index_of_nearest_xy(lon, lat, aeri_lon, aeri_lat)
    size = int(size)
    temp =  d.variables['Temperature'][0,:,idy-size:idy+size,idx-size:idx+size]
    rh = d.variables['Relative_humidity'][0,:,idy-size:idy+size,idx-size:idx+size]
    hght = d.variables['Geopotential_height'][0,:,idy-size:idy+size,idx-size:idx+size]
    
    sfc_pres = d.variables['Pressure_surface'][0,idy-size:idy+size,idx-size:idx+size]
    sfc_hght = d.variables['Geopotential_height_surface'][0,idy-size:idy+size,idx-size:idx+size]
    sfc_temp = d.variables['Temperature_height_above_ground'][0,0,idy-size:idy+size,idx-size:idx+size]
    sfc_rh = d.variables['Relative_humidity_height_above_ground'][0,0,idy-size:idy+size,idx-size:idx+size]
    
    d.close()
    
    mxrs = []
    temps = []
    press = []
    rhs = []

    for index, x in np.ndenumerate(sfc_hght):
        idx_aboveground = np.where(sfc_hght[index] < hght[:,index[0], index[1]])[0]

        new_hght = (np.hstack((sfc_hght[index]+2, hght[idx_aboveground, index[0], index[1]])) - sfc_hght[index])/1000.
        new_temp = np.hstack((sfc_temp[index], temp[idx_aboveground, index[0], index[1]]))
        new_rh = np.hstack((sfc_rh[index], rh[idx_aboveground, index[0], index[1]]))
        new_pres = np.hstack((sfc_pres[index], pres[idx_aboveground]))
                
        new_q = rh2q(new_temp, new_pres, new_rh/100.)*1000.
        
        oe_mxr = np.interp(aerioe_hght, new_hght, new_q)
        oe_temp = np.interp(aerioe_hght, new_hght, new_temp)
        oe_pres = np.interp(aerioe_hght, new_hght, new_pres)
        oe_rh = np.interp(aerioe_hght, new_hght, new_rh)
        mxrs.append(oe_mxr)
        temps.append(oe_temp)
        press.append(new_pres)
        rhs.append(oe_rh)

    temps = np.asarray(temps)
    mxrs = np.asarray(mxrs)
    rhs = np.asarray(rhs)

    return temps, mxrs, press, type, path
"""

# This function is used to find the nearest point on the grid to the AERI location
def find_index_of_nearest_xy(y_array, x_array, y_point, x_point):
    distance = (y_array-y_point)**2 + (x_array-x_point)**2
    idy,idx = np.where(distance==distance.min())
    return idy[0],idx[0]

# This is used to convert between the relative humidity grid that is stored in the
# model grids to mixing ratio (which is what we need for the prior)
def rh2q(temp, pres, rh):
    Rv = 461.
    L = 2.453 * 10**6
    es = 6.11 * np.exp((L/Rv)*((1./(273.15)) - (1./temp)))
    e = rh * es
    q = (0.622*e) / ((pres/100.) - e)
    return q


# This the code that gets called by run_prior_gen.py when we want to make a prior out of the
# ARM RUC/RAP MODEL DATA.
def getARMModelPrior(model_data_path, climo_prior, yyyymmdd, size, hh, hh_delta, aeri_lon, aeri_lat):
    print "This prior is spatially centered at: " + str(aeri_lat) + ',' + str(aeri_lon)

    climo = Dataset(climo_prior)
    aerioe_hght = climo.variables['height'][:]

    dt = datetime.strptime(yyyymmdd + hh, '%Y%m%d%H') - timedelta(seconds=int(hh_delta)*60*60)
    end_dt = datetime.strptime(yyyymmdd + hh, '%Y%m%d%H') + timedelta(seconds=int(hh_delta)*60*60)
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
        temp, mxr, pres, type, link = getARMProfiles(model_data_path, yyyymmdd, hour, aeri_lon, aeri_lat, size, aerioe_hght)
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
def getRealtimePrior(climo_prior, yyyymmdd, size, hh, hh_delta, aeri_lon, aeri_lat):
    #Load in latitude and longitude points for the RAP/RUC 13 km grid.

    print "Will be pulling observations from the nearest point to: " + str(aeri_lat) + ', ' + str(aeri_lon)
    
    #Get the nearest grid point to the AERI location.
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


# This is the code that gets called by run_prior_gen.py when we want to make a prior from
# ONLINE RUC/RAP data.
def getOnlineModelPrior(climo_prior, yyyymmdd, size, hh, hh_delta, aeri_lon, aeri_lat):
    #Load in latitude and longitude points for the RAP/RUC 13 km grid.

    print "This prior is centered at: " + str(aeri_lat) + ',' + str(aeri_lon)
    
    #Get the nearest grid point to the AERI location.

    climo = Dataset(climo_prior)
    aerioe_hght = climo.variables['height'][:]

    dt = datetime.strptime(yyyymmdd + hh, '%Y%m%d%H') - timedelta(seconds=int(hh_delta)*60*60)

    all_temps = None
    all_mxrs = None
    all_pres = None
    print "Gathering profiles within a " + str(2*size) + "x" + str(2*size) + " grid."
    types = []
    paths = []
    for i in range(0, int(hh_delta)*2,1):
        dt_string = datetime.strftime(dt, '%Y%m%d%H')
        forecast_hour = (3-len(str(i)))*'0' + str(i) 
        print "Gathering profiles for model run on " + dt_string + " F" + forecast_hour
        # This only searches the NOMADS THREDDS Server
        temp, mxr, pres, type, link = getHourlyProfiles(dt_string, forecast_hour, aeri_lat, aeri_lon, size, aerioe_hght)
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
    sfc_mutiplier = 3
    yyyymmdd = '20130531'
    hour = '12'
    print "Testing RUC time frame"
    mean, cov, climo, type, path, yyyymmdd, hh, n = getModelPrior('prior_data/Xa_Sa_datafile.55_levels.month_05.cdf', yyyymmdd, 10, hour, 3) 
    fig = figure(figsize=(12,6))
    subplot(121)
    title("W/O Inflation " + yyyymmdd + ' ' + hour)
    contourf(cov, np.arange(-14,15,1), cmap=get_cmap('RdBu'), extend='both')
    axvline(x=55, c='k')
    axhline(y=55, c='k')
    
    top = 3
    d = Dataset('prior_data/Xa_Sa_datafile.55_levels.month_05.cdf')
    height = d.variables['height'][:]
    d.close()
    profile_type = 2
    new_cov = inflatePrior(cov, sfc_mutiplier, top, height, profile_type)
    subplot(122)
    title("Inflation " + yyyymmdd + ' ' + hour)
    im = contourf(new_cov, np.arange(-14,15,1), cmap=get_cmap('RdBu'), extend='both')
    axvline(x=55, c='k')
    axhline(y=55, c='k')
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)
    show()
    stop
    print "Testing RUC time frame"
    getModelPrior('prior_data/Xa_Sa_datafile.55_levels.month_08.cdf', "20110103", 10, "12", 3)

if __name__ == '__main__':
    aeri_lon = -97
    aeri_lat = 35
    print getMotherlodeProfiles('2015022412', 6, aeri_lat, aeri_lon, 10, [10,20, 30])
    #main()

