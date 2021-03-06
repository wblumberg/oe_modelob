##Model Observation Script
Author: Greg Blumberg OU/CIMMS

Recent changes to AERIoe by Dave Turner have included the option to have a thermodynamic profile be a "Y-vector" observation in the retrieval.  This facilitates the addition of radiosonde/Raman Lidar/DIAL data as inputs into the AERIoe retrieval.  As a consequence, we have decided to remove the AERI model prior concept from our processing and instead use a model thermodynamic profile as another "observation" into the AERIoe retrieval.  This enables a better retrieval in the mid-to-upper atmosphere.  This script will perform the conversions from the model grids into the model-profile netCDF files that AERIoe can now accept as inputs. 

As of 28 Oct 2016, both the ARM and Motherlode data sources also support the download and storing of the model wind speed and direction profile.  The lowest data point for these profiles should be 10 m AGL, and any missing data is stored in the output netCDF file as the masked or missing field.  NOTE: The model U and V components are linearly interpolated to a much finer vertical grid before calculating the wind speed and direction.

This script takes in a set of either real-time or archived model data (e.g. ARM-RUC grids or realtime RAP/HRRR data or NCDC archived RUC data or realtime GFS) and a lat/lon point.  It then pulls the temperature and water vapor mixing ratio profile for each time at that point to develop a netCDF file.  The output file contains the time/height data for temperature and water vapor mixing ratio, pressure, height, temp_std, and wvmr_std. It also contains metadata about how the file was made.

      This code currently supports making these files from:
      - realtime RAP profiles
      - ARM-formatted RUC/RAP files

      Requires:
            netCDF4-python (conda install netcdf4)

      Arguments:
          [1] YYYYMMDD - date of AERI file to run (ex: 20130531)
          [2] VIP file - path to config file (ex: nwc_vip.txt)
          [3] BHOUR - beginning hour of profiles (ex: 5)
          [4] EHOUR - ending hour of profiles (ex: 5)

       Output:
           Files of the name format of: RRmodelsoundings.YYYYMMDD.HH.LAT.LON.cdf

      Example command line arguments:
            $ python run_ModelOB.py 20150306 vip.txt 00 00 (creates a file with only the 00 UTC model observation)
            $ python run_ModelOB.py 20030508 vip.txt 00 23 (creates a file for the day of 2003-05-08)

#### VIP File Variables:

The VIP file is similar to the AERIoe format in order to reduce the number of files needed to run AERIoe.
Here are the variables included in the VIP file, along with a description of them:
   
      Variable: data_source
      Required?: YES
      Options: 1 - NOAA NOMADS RAP/RUC analysis archive
               2 - ARM-formatted RAP/RUC hourly analysis files
               3 - 13 km RAP Motherlode UCAR data server (forecast and analysis data)
      Description:
      These three options can be used to generate files for different AERI data.
      
         1 - contains analysis data for the CONUS back to 2010, however data outages are sporatic.  
         More recent data may be less sporadic. This data source is recommended when you need model 
         observations outside of the SGP, however it may work for the 2011-2013 NWC AERI deployment.
         
         Browse this datastream at this website:
         http://nomads.ncdc.noaa.gov/thredds/catalog/rap130/
      
         2 - Contains analysis data for the SGP back to the 90s.  You'll have to order this data 
         from the ARM archive however data for this does exist on the ARM cluster and the RAID.  The
         data files need to be hourly and of the format:
      
                  sgpallruc60X1.c0        05/08/1996  to  04/23/1998  (not tested)
                  sgpallruc40isobX1.c1    04/20/1998  to  04/16/2002  (tested)
                  sgpsynruc20isobX1.c1    04/18/2002  to  05/01/2012  (tested)
                  sgpsynrap20plevX1.c1    05/01/2012  to  current     (tested)
      
              The ARM archive also supports data on the native grid (hybrid levels).
              Support for this data type has not been implemented in this program, however it
              may happen in the future since the hybrid levels provide a higher resolution data grid than
              these isobaric files.
              
              More information about this data type can be found at:
              http://www.arm.gov/xdc/xds/ruc
      
              Selection of this data source option requires implementing the "arm_model_dir" VIP variable.
      
          3 - Uses analysis and forecast data from the 13 km RAP data located on the UCAR Motherlode server.
          This data source semi-realtime data from the last 2 weeks and is updated regularly.  When using
          this data source, the program will use forecast data if analysis data does not exist for the
          time frame specified.  It is recommended that this data source is used for real-time generation
          of RAP Model Observation files.  
          
          Browse this datastream at:
          http://thredds.ucar.edu/thredds/catalog/grib/NCEP/RAP/CONUS_13km/catalog.html
      
      
      Variable: aeri_lat
      Required?: YES
      Options: None (float; decimal degrees)
      Description:  This variable is the latitude value of where you'd like the 
      profile center at (e.g. the AERI location)
      
      Variable: aeri_lon
      Required?: YES
      Options: None (float; decimal degrees)
      Description:  This variable is the longitude value of where you'd like the 
      profile center at       (e.g. the AERI location)
      
      Variable: arm_model_dir
      Required?: only if data_source=2
      Options: None (string type)
      Description:  The location of the ARM-formatted RAP or RUC data.
      
      Variable: spatial_mesh_size
      Required?: YES
      Options: None (integer)
      Description:  Half of the number of grid points along the permeter of the spatial mesh 
      (e.g. a value of 5 is a 10x10 mesh centered near aeri_lat, aeri_lon.)
      
      Variable: temporal_mesh_size
      Required?: YES
      Options: None (integer)
      Description:  The number of hours of the temporal data window to include in the std calculation 
      (e.g. a value of 2 means the STD calculation will use 5 hours of data.)
      
      Variable: model_ob_output_dir
      Required?: YES
      Options: None (string)
      Description:  The directory the model observation files should be written to.
      
#### Troubleshooting

1. If your selected date/time range has no data, the program should break with an error message saying that there was no data.  It is helpful to check the data source websites (see above under VIP File Variables) to see if the data actually exists on the server.  If it isn't, you're SOL until I can make some modifications to include other data types for the program to use.
