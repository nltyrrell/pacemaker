import numpy as np
import numpy.ma as ma
import iris as iris
import cartopy.crs as ccrs
import scipy.stats as stats
import sys

Rs_da=287.05          # Specific gas const for dry air, J kg^{-1} K^{-1}
Cp_da=1004.6          # Specific heat at constant pressure for dry air

def theta(temp_cube,RH_cube,tempC_cube,pres=False,pref=1000.):
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
    tempC_data = tempC_cube.data
    RH_data = RH_cube.data
    pres_hPa = np.ones(temp_data.shape)
    pres = np.expand_dims(np.expand_dims(pres,axis=1),axis=1) 
    pres_hPa = pres*pres_hPa

    theta = temp_data*(pref/pres)**(Rs_da/Cp_da)
    print pres
    print pref

    r_L = 0
    # Temp where above/below 0C
    TC_ice = np.ma.masked_greater_equal(tempC_data,0.0)
    TC_liq = np.ma.masked_less(tempC_data,0.0)

    # Calculate sat vap pres, e_s[hPa]
    e_s_liq_mask = 6.112*np.exp(17.67*TC_liq/(TC_liq+243.12))
    e_s_ice_mask = 6.112*np.exp(22.46*TC_ice/(TC_ice+272.62))
    e_s_ice = e_s_ice_mask.data*~e_s_ice_mask.mask
    e_s_liq = e_s_liq_mask.data*~e_s_liq_mask.mask
    e_s = e_s_liq + e_s_ice 

    # Calc Sat MixRatio r_s[kg/kg]
    r_s = 0.622*(e_s/(pres_hPa-e_s))

    # Calculate mix ratio [kg/kg]
    r = r_s*RH_data/100

    theta_v = theta * (1 + 0.61*r - r_L)

    return theta, theta_v, r, r_s


ncfile_path = '/home/nicholat/project/pacemaker/ncfiles/'
# temp_plv = iris.load_cube(ncfile_path + 'temp.plv.4ysl.m48.nc')
temp_abs = iris.load_cube(ncfile_path + 'temp.plv.4ysl.m48.abs.nc')
temp_abs.coord('p').standard_name = 'air_pressure'
tempC = temp_abs.copy()
tempC.convert_units('celsius')

RH_abs = iris.load_cube(ncfile_path + 'rhum.plv.4ysl.m48.abs.nc')
# temp_plv.coord('p').standard_name = 'air_pressure'

theta, theta_v, r, r_s = theta(temp_abs,RH_abs,tempC)

theta_cube = temp_abs.copy()
thetav_cube = temp_abs.copy()

theta_cube.data = theta
theta_cube.standard_name = 'air_potential_temperature'
# theta_cube.units = 'K'
theta_cube.attributes['title'] = 'Theta'
theta_cube.attributes['name'] = 'theta'
# theta_cube.remove_coord('surface')
# theta_cube.remove_coord('time')
theta_mean = theta_cube.collapsed('time',iris.analysis.MEAN)
theta_anom = theta_cube - theta_mean
iris.save(theta_anom,ncfile_path+'theta.plv.4ysl.m48.nc')

thetav_cube.data = theta_v
thetav_cube.standard_name = 'virtual_temperature'
# thetav_cube.units = 'K'
thetav_cube.attributes['title'] = 'ThetaV'
thetav_cube.attributes['name'] = 'thetav'
# thetav_cube.remove_coord('surface')
# thetav_cube.remove_coord('time')
thetav_mean = thetav_cube.collapsed('time',iris.analysis.MEAN)
thetav_anom = thetav_cube - theta_mean
iris.save(thetav_anom,ncfile_path+'thetav.plv.4ysl.m48.nc')













