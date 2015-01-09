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
import prettyplotlib as ppl
import pickle

"""
Test the enso mask, for use in oc_ln_varresponse.py
"""

ncfile_path = '/home/nicholat/project/pacemaker/ncfiles/'

temp_plv = iris.load_cube(ncfile_path + 'temp.plv.4ysl.m48.nc')
# temp_reg = regions(temp_plv)
temp_sfc = iris.load_cube(ncfile_path + 'temp.sfc.4ysl.m48.nc')
# stemp_reg = regions(temp_sfc)
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
#     print 'mean sfc '+str(cube_regmean[0].data)
#     print 'mean p=10 '+str(cube_regmean[10].data)
    return cube_reg, cube_regmean


# def regions(cube,clim=False,maxf=True):
clim=False;maxf=True
cube = temp_sfc
mons = 4
lag = 1
max_i = 35 + lag
max_f = max_i + mons
min_i = 11 + lag
min_f = min_i + mons
if maxf:
    forc_i = max_i
    forc_f = max_f
else:
    forc_i = min_i
    forc_f = min_f


if clim:
    cube_max = cube.collapsed('time',iris.analysis.MEAN)
else:
    cube_max = cube[forc_i:forc_f,:,::].collapsed('time',iris.analysis.MEAN)
# cube_anom = cube-cube_mean

cube_max.coord('latitude').guess_bounds()
cube_max.coord('longitude').guess_bounds()

print cube_max.shape
# ---------- define sst/tland----------
lsmask = iris.load_cube(ncfile_path + 'lsmask.nc')[0,0,::]
landmask = ~(ma.make_mask(lsmask.data.copy()) + np.zeros(cube_max.shape)).astype(bool) # mask sea, show land

# mask out the enso region: masked = 1, not masked = 0
# mask from x = 43:73, y = 32:40
enso_mask = np.zeros(cube_max.shape)
enso_mask[0,30:42,43:75] = 1
seamask = (ma.make_mask(lsmask.data.copy()) + np.zeros(cube_max.shape) + enso_mask).astype(bool) # mask land, show sea

Cocean = cube_max.copy()
Cland = cube_max.copy()
Cocean.data = ma.array(Cocean.data, mask=seamask)
Cland.data = ma.array(Cland.data, mask=landmask)
# --------------
land_cube = Cland
ocean_cube = Cocean

oclat10, oclat10_mean = regmean(ocean_cube,loni=0,lonf=360,lati=-10,latf=10)
oclatp20, oclatp20_mean = regmean(ocean_cube,loni=0,lonf=360,lati=10,latf=20)
oclatm20, oclatm20_mean = regmean(ocean_cube,loni=0,lonf=360,lati=-20,latf=-10)
oclatp30, oclatp30_mean = regmean(ocean_cube,loni=0,lonf=360,lati=20,latf=30)
oclatm30, oclatm30_mean = regmean(ocean_cube,loni=0,lonf=360,lati=-30,latf=-20)

lat10, lat10_mean = regmean(land_cube,loni=0,lonf=360,lati=-10,latf=10)
latp20, latp20_mean = regmean(land_cube,loni=0,lonf=360,lati=10,latf=20)
latm20, latm20_mean = regmean(land_cube,loni=0,lonf=360,lati=-20,latf=-10)
latp30, latp30_mean = regmean(land_cube,loni=0,lonf=360,lati=20,latf=30)
latm30, latm30_mean = regmean(land_cube,loni=0,lonf=360,lati=-30,latf=-20)

