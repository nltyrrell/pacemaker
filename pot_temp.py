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


ncfile_path = '/home/nicholat/project/pacemaker/ncfiles/'
temp_plv = iris.load_cube(ncfile_path + 'temp.m48.plv.4ysl.nc')
temp_plv.coord('p').standard_name = 'air_pressure'

theta_plv = theta(temp_plv)
theta_cube = temp_plv.copy()

theta_cube.data = theta_plv
theta_cube.long_name = 'Theta'
theta_cube.units = 'K'
theta_cube.attributes['title'] = 'Theta'
theta_cube.attributes['name'] = 'theta'
# theta_cube.remove_coord('surface')
# theta_cube.remove_coord('time')
iris.save(theta_cube,ncfile_path+'theta.plv.4ysl.m48.nc')
sys.exit('exitttt')
def theta(temp_cube,RH_cube,pres=False,pref=100000.):
    """ Calculate the Virtual Potential Temperature
    theta_v = theta * (1 + 0.61*r - r_L)
    Cloudy/sat air use r = r_s
    Unsat air use r_L = 0
    to find r
    r = r_s*RH/100
    r_s = 6.112*exp(22.46*tempC/(tempC+272.62))*100.

    Input:
    temp_cube [K]
    Output: Theta [K]
    """
    pres = temp_cube.coord('air_pressure').points
    temp_data = temp_cube.data
    new_pres = np.ones(temp_data.shape)
    pres = np.expand_dims(np.expand_dims(pres,axis=1),axis=1) 
    new_pres = pres*new_pres

    theta = temp_data*(pref/pres)**(Rs_da/Cp_da)

    tempC = temp_data - 273.15

    r_s = 6.112*exp(22.46*tempC/(tempC+272.62))*100.

    theta_v = theta * (1 + 0.61*r - r_L)

    return theta, theta_v













