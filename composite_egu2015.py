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

ncfile_path = '/home/nicholat/project/pacemaker/ncfiles/'
temp = iris.load_cube(ncfile_path + 'temp.sfc.4ysl.m48.nc')
# pres = iris.load_cube(ncfile_path + 'pres.sfc.4ysl.m48.nc')

# import the ACCESS data using iris
temp_plv = iris.load_cube(ncfile_path + 'temp.plv.4ysl.m48.nc')
temp_plv.coord('p').standard_name = 'air_pressure'

# Define regions
temp_plv.coord('latitude').guess_bounds()
temp_plv.coord('longitude').guess_bounds()
# pres.coord('latitude').guess_bounds()
# pres.coord('longitude').guess_bounds()

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

mons = 3
lag = 3
max_i = 35 + lag
max_f = max_i + mons
min_i = 11 + lag
min_f = min_i + mons

temp_maxT_mean  = temp[max_i:max_f,0,::].collapsed('time',iris.analysis.MEAN)
Tmax            = temp_maxT_mean
temp_minT_mean  = temp[min_i:min_f,0,::].collapsed('time',iris.analysis.MEAN)
Tmin            = temp_minT_mean

temp_300            = temp_plv.extract(iris.Constraint(air_pressure=300))
temp_300_maxT_mean  = temp_300[max_i:max_f,::].collapsed('time',iris.analysis.MEAN)
T300max             = temp_300_maxT_mean
temp_300_minT_mean  = temp_300[min_i:min_f,::].collapsed('time',iris.analysis.MEAN)
T300min             = temp_300_minT_mean

plt.close('all')
fig1 = plt.figure(1)
ax = plt.subplot(projection=ccrs.Robinson())
mesh = iplt.pcolormesh(Tmax,vmin=-1,vmax=1,cmap=mc.jetwhite())
ax.coastlines()
plt.colorbar(mesh, orientation='horizontal',extend=
        'both',label='K')
plt.title('T response to max positive forcing')
plt.show()
plt.savefig('figures/comp_Tmax_sfc.png')

fig2 = plt.figure(2)
ax = plt.subplot(projection=ccrs.Robinson())
mesh = iplt.pcolormesh(T300max,vmin=-1,vmax=1,cmap=mc.jetwhite())
ax.coastlines()
plt.title('T response to max positive forcing, 300hPa')
plt.colorbar(mesh, orientation='horizontal',extend=
        'both',label='K')
plt.show()
plt.savefig('figures/comp_Tmax_300.png')


# plt.figure(3)
# qplt.pcmeshclf(T300max - T300min,vmin=-1,vmax=1,cmap=mc.jetwhite())
# plt.savefig('figures/reg_T_sfc_RH_sfc.png')









