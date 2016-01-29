import numpy as np

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