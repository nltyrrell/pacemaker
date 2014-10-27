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


# Load in Tsfc, Ttropos (or T300hPa?)
# calculate the regression between the timeseries
# Tsfc = alpha * Ttropo + epsilon

ncfile_path = '/home/nicholat/project/mit_tcm/access_runs/ncfiles/'
temp = iris.load_cube(ncfile_path + 'temp.m48.sfc.4ysl.nc')

# import the ACCESS data using iris
temp_plv = iris.load_cube(ncfile_path + 'temp.m48.plv.4ysl.nc')
temp_plv.coord('p').standard_name = 'air_pressure'

# Define regions
temp_plv.coord('latitude').guess_bounds()
temp_plv.coord('longitude').guess_bounds()

# ---------- define sst/tland----------
lsmask = iris.load_cube(ncfile_path + 'lsmask.nc')[0,0,::]
landmask = ~(ma.make_mask(lsmask.data.copy()) + np.zeros(temp_plv.shape)).astype(bool) # mask sea, show land
seamask = (ma.make_mask(lsmask.data.copy()) + np.zeros(temp_plv.shape)).astype(bool) # mask land, show sea

# sys.exit('exiiiiiiitt')
Tocean = temp_plv.copy()
Tland = temp_plv.copy()
Tocean.data = ma.array(Tocean.data, mask=seamask)
Tland.data = ma.array(Tland.data, mask=landmask)
# --------------
tropics = iris.Constraint(latitude = lambda v: -30 <= v <= 30)
Tocean_trop = Tocean.extract(tropics)
grid_areas_trop = iris.analysis.cartography.area_weights(Tocean_trop)

temp_maxT_mean  = temp[35:41,0,::].collapsed('t',iris.analysis.MEAN)
Tmax            = temp_maxT_mean
temp_minT_mean  = temp[11:17,0,::].collapsed('t',iris.analysis.MEAN)
Tmin            = temp_minT_mean

temp_300            = temp_plv.extract(iris.Constraint(air_pressure=300))
temp_300_maxT_mean  = temp_300[35:41,::].collapsed('t',iris.analysis.MEAN)
T300max             = temp_300_maxT_mean
temp_300_minT_mean  = temp_300[11:17,::].collapsed('t',iris.analysis.MEAN)
T300min             = temp_300_minT_mean


reg_cube = copy_cube.copy()
reg_cube.data[:] = linreg_map
reg_cube.long_name = 'Lin Regression '+ name1 +' '+ name2
reg_cube.units = 'no_unit'
reg_cube.attributes['title'] = 'Lin Regression '+ name1 +' '+ name2
reg_cube.attributes['name'] = 'reg'
reg_cube.remove_coord('surface')
reg_cube.remove_coord('time')
iris.save(reg_cube,ncfile_path+'lreg.4ysl.'+name1+'.'+name2+'.nc')







