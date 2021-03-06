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

mons = 3
lag = 3
max_i = 35 + lag
max_f = max_i + mons
min_i = 11 + lag
min_f = min_i + mons

temp_maxT_mean  = temp[max_i:max_f,0,::].collapsed('t',iris.analysis.MEAN)
Tmax            = temp_maxT_mean
temp_minT_mean  = temp[min_i:min_f,0,::].collapsed('t',iris.analysis.MEAN)
Tmin            = temp_minT_mean

temp_300            = temp_plv.extract(iris.Constraint(air_pressure=300))
temp_300_maxT_mean  = temp_300[max_i:max_f,::].collapsed('t',iris.analysis.MEAN)
T300max             = temp_300_maxT_mean
temp_300_minT_mean  = temp_300[min_i:min_f,::].collapsed('t',iris.analysis.MEAN)
T300min             = temp_300_minT_mean

temp_700            = temp_plv.extract(iris.Constraint(air_pressure=700))
temp_700_maxT_mean  = temp_700[max_i:max_f,::].collapsed('t',iris.analysis.MEAN)
T700max             = temp_700_maxT_mean
temp_700_minT_mean  = temp_700[min_i:min_f,::].collapsed('t',iris.analysis.MEAN)
T700min             = temp_700_minT_mean

plt.close('all')
plt.ion()
# plt.figure(1)
qplt.pcmeshclf(Tmax-T700max,vmin=-0.5,vmax=0.5,cmap=mc.jetwhite())
plt.title('T response to max positive forcing')
# plt.savefig('figures/comp_Tmax_sfc.png')

# plt.figure(2)
# qplt.pcmeshclf(Tmin,vmin=-1,vmax=1,cmap=mc.jetwhite_r())
# plt.title('T response to max negative forcing')
# plt.savefig('figures/comp_Tmin_sfc.png')
# 
# plt.figure(3)
# qplt.pcmeshclf(T300max,vmin=-1,vmax=1,cmap=mc.jetwhite())
# plt.title('T response to max positive forcing, 300hPa')
# plt.savefig('figures/comp_Tmax_300.png')
# 
# plt.figure(4)
# qplt.pcmeshclf(T300min,vmin=-1,vmax=1,cmap=mc.jetwhite_r())
# plt.title('T response to max negative forcing, 300hPa')
# plt.savefig('figures/comp_Tmin_300.png')
# 
# plt.figure(5)
# qplt.pcmeshclf(T850max,vmin=-1,vmax=1,cmap=mc.jetwhite())
# plt.title('T response to max positive forcing, 850hPa')
# plt.savefig('figures/comp_Tmax_850.png')
# 
# plt.figure(6)
# qplt.pcmeshclf(T850min,vmin=-1,vmax=1,cmap=mc.jetwhite_r())
# plt.title('T response to max negative forcing, 850hPa')
# plt.savefig('figures/comp_Tmin_850.png')

# plt.figure(3)
# qplt.pcmeshclf(T300max - T300min,vmin=-1,vmax=1,cmap=mc.jetwhite())
# plt.savefig('figures/reg_T_sfc_RH_sfc.png')









