import numpy as np
import sys

# This function performs a unit check on the profiles and makes sure they are correct.
def unitCheck(prof_data):
    print("Performing unit tests on the data collected by the program...")
    # Check the pressure value
    if np.ma.max(prof_data['pressure']) > 10000:
        # Convert from Pascals to mb
        print("\tThe values in the pressure array are greater than 10000, let's convert to mb from Pascals...")
        prof_data['pressure'] = prof_data['pressure'] / 100. 
#    if len(np.ma.where(prof_data['wvmr'] < 0.01)[1]) / float(len(prof_data['wvmr'][0])) > .7:
#        print "\tThere might be a unit error for the water vapor mixing ratio profile...the values are too small...converting..."
#        prof_data['wvmr'] = prof_data['wvmr'] * 1000.
#        prof_data['wvmr_sigma'] = prof_data['wvmr_sigma'] * 1000.
    #if np.ma.max(prof_data['temperature']) > 100:
    #    print "\tA temperature value is above 100 C...this is probably a unit error and may be affecting the WVMR calculation..."
    #    print "\tThey've changed something in the Matrix...I mean data server...very likely need to rewrite this script..."
    #    print "FAILED."
    #    sys.exit()
    return prof_data

# This function is used to find the nearest point on the grid to the AERI location
def find_index_of_nearest_xy(y_array, x_array, y_point, x_point):
    distance = (y_array-y_point)**2 + (x_array-x_point)**2
    idy,idx = np.where(distance==distance.min())
    return idy[0],idx[0]

# This is used to convert between the relative humidity grid that is stored in the
# model grids to mixing ratio (which is what we need for the prior)
def rh2q(temp, pres, rh):
    """
        Inputs
        ------
        temp [K]
        pres [Pa]
        rh [fraction]
    """
    Rv = 461.
    L = 2.453 * 10**6
    es = 6.11 * np.exp((L/Rv)*((1./(273.15)) - (1./temp)))
    e = rh * es
    q = (0.622*e) / ((pres/100.) - e)
    return q

# This function is used to convert from sigma/hybrid coordinates to pressure coordinates.
def get3Dpres(hbcoefa, hbcoefb, P0, sfc_pres):
    # hbcoefa(k) ( vertical coord)
    # hbcoefb(k) ( vertical coord)
    # P0 - float
    # sfc_pres(i,j) (2-d grid)
    # Returns P(i,j,k)
    pres_arr = np.empty(sfc_pres.shape)
    for i in range(len(pres_arr)):
        pres_arr[i] = (hbcoefa[i]*P0) + hbcoefb[i] * sfc_pres[i]
    return pres_arr

def hypsometric_z(pres_arr, temp_arr, q_arr, sfc_z):
    new_z = np.empty(pres_arr.shape)
    new_z[0,:,:] = sfc_z
    for i in np.ndenumerate(sfc_z):
        for j in range(1, len(pres_arr)):
            x_idx = i[0][0]
            y_idx = i[0][1]
            
            mean_t = (temp_arr[j, x_idx, y_idx] + temp_arr[j-1, x_idx, y_idx])/2.
            z_1 = new_z[j-1,x_idx,y_idx]
            p_1 = pres_arr[j-1,x_idx,y_idx]
            p_2 = pres_arr[j,x_idx,y_idx]
            new_z[j,x_idx,y_idx] = z_1 + (np.log(p_1/p_2) * (287. * mean_t) / (9.81))

    return new_z




