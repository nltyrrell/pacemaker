import numpy as np
import numpy.ma as ma
import iris as iris
import iris.plot as iplt
import iris.quickplot as qplt
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import mycmaps as mc
import scipy.stats as stats
import sys
import troposave as ta

"""
Define many regions, based on surface temperature response (or moisture?)
get the area mean value in these regions for various variables
plot... in some useful mannner tbc
"""

def regmean(cube,loni,lonf,lati,latf):
    """ Define a region and get the area weighted mean
    Input:  cube, lon_i, lon_f, lat_i, lat_f
    Output: cube_reg, cube_regmean
    """
    region = iris.Constraint(longitude=lambda l: (loni <= l <= lonf), latitude = lambda l: (lati <= l <= latf))
    cube_reg = cube.extract(region)
    grid_areas = iris.analysis.cartography.area_weights(cube_reg)
    cube_regmean = cube_reg.collapsed(['latitude', 'longitude'],
                           iris.analysis.MEAN,
                           weights=grid_areas)
    print 'mean sfc '+str(cube_regmean[0].data)
#     print 'mean p=10 '+str(cube_regmean[10].data)
    return cube_reg, cube_regmean


def regions(cube,clim=False):
    if clim:
        cube_max = cube.collapsed('time',iris.analysis.MEAN)
    else:
        cube_max = cube[35:42,:,::].collapsed('time',iris.analysis.MEAN)
    # cube_anom = cube-cube_mean

    cube_max.coord('latitude').guess_bounds()
    cube_max.coord('longitude').guess_bounds()

    # ---------- define sst/tland----------
    lsmask = iris.load_cube(ncfile_path + 'lsmask.nc')[0,0,::]
    landmask = ~(ma.make_mask(lsmask.data.copy()) + np.zeros(cube_max.shape)).astype(bool) # mask sea, show land
    seamask = (ma.make_mask(lsmask.data.copy()) + np.zeros(cube_max.shape)).astype(bool) # mask land, show sea

    Cocean = cube_max.copy()
    Cland = cube_max.copy()
    Cocean.data = ma.array(Cocean.data, mask=seamask)
    Cland.data = ma.array(Cland.data, mask=landmask)
    # --------------
    land_cube = Cland

    India, India_mean = regmean(land_cube,loni=60,lonf=75,lati=0,latf=30)
    MC, MC_mean = regmean(land_cube,loni=90,lonf=140,lati=-10,latf=10)
    TropSthAm, TropSthAm_mean = regmean(land_cube,loni=290,lonf=315,lati=-23,latf=0)
    SthSthAm, SthSthAm_mean = regmean(land_cube,loni=270,lonf=315,lati=-60,latf=-24)
    NthWestAfr, NthWestAfr_mean = regmean(land_cube,loni=-15,lonf=15,lati=10,latf=30)
    NthEastAfr, NthEastAfr_mean = regmean(land_cube,loni=15,lonf=50,lati=10,latf=30)
    TropAfr, TropAfr_mean = regmean(land_cube,loni=12,lonf=40,lati=-15,latf=5)
    Aus, Aus_mean = regmean(land_cube,loni=120,lonf=140,lati=-30,latf=-17)

    return [India_mean , MC_mean , TropSthAm_mean , SthSthAm_mean , NthWestAfr_mean , NthEastAfr_mean , TropAfr_mean , Aus_mean]
    
ncfile_path = '/home/nicholat/project/pacemaker/ncfiles/'

temp_plv = iris.load_cube(ncfile_path + 'temp.plv.4ysl.m48.nc')
temp_reg = regions(temp_plv)
rhum_plv = iris.load_cube(ncfile_path + 'rhum.plv.4ysl.m48.nc')
rhum_reg = regions(rhum_plv)
smc = iris.load_cube(ncfile_path + 'smc.sfc.4ysl.m48.nc')
smc_clim_reg = regions(smc,clim=True)
smc_reg = regions(smc)

# plt.clf()
# plt.plot(lnochold,pressure,'b')#color=wetdry(k,nloops))
# plt.legend(('NthAfr','Central Africa','Sth America','Australia','Tropical Ocean'),loc=3, )
# plt.title('Regression with surface temp and temp profile.') 
# plt.xlim(-0.5,3.0)
# plt.ylabel('z [hPa]')
# plt.xlabel('Regresion coeff.')
# plt.gca().invert_yaxis()
# plt.savefig('./figures/regprof_cont.png')




