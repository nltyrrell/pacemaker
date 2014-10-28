import numpy as np
import numpy.ma as ma
import iris as iris
import cartopy.crs as ccrs
import scipy.stats as stats
import sys

Rs_da=287.05          # Specific gas const for dry air, J kg^{-1} K^{-1}
Cp_da=1004.6          # Specific heat at constant pressure for dry air

def theta(temp_cube,pres=False,pref=100000.):
    """ Calcualte Potential Temperature
    Input:
    temp_cube [K]
    Output: Theta [K]
    """
    p0 = iris.coords.AuxCoord(1000.0,
                              long_name='reference_pressure',
                              units='hPa')
    pres = temp_cube.coord('air_pressure').points
    temp_data = temp_cube.data
    new_pres = np.ones(temp_data.shape)
    pres = np.expand_dims(np.expand_dims(pres,axis=1),axis=1) 
    new_pres = pres*new_pres

    pot_temp = temp_data*(pref/pres)**(Rs_da/Cp_da)

    return pot_temp
