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

pressure = temp_plv.coord('air_pressure').points

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
TO_trop_mean = Tocean_trop.collapsed(['latitude','longitude'],
                iris.analysis.MEAN,weights=grid_areas_trop)

Tocean_forc = Tocean.extract(iris.Constraint(latitude = lambda v: -20 <= v <= 20,
                longitude = lambda l: 140 <= l <= 300)).collapsed(['latitude','longitude'],
                iris.analysis.MEAN,weights=iris.analysis.cartography.area_weights(Tocean_forc))

trop_remote = 
Tocean_rem = Tocean.extract(iris.Constraint(latitude = lambda v: -30 <= v <= 30,
                longitude = lambda l: (300 <= l <= 360) | (0 <= l <= 130) )).collapsed(['latitude','longitude'],
                iris.analysis.MEAN,weights=iris.analysis.cartography.area_weights(Tocean_rem))

Tland_trop = Tland.extract(tropics)
grid_areas_trop = iris.analysis.cartography.area_weights(Tland_trop)
TL_trop_mean = Tland_trop.collapsed(['latitude','longitude'],
                iris.analysis.MEAN,weights=grid_areas_trop)
mons = 4
lag = 3
max_i = 35 + lag
max_f = max_i + mons
min_i = 11 + lag
min_f = min_i + mons
TOf_max  = TO_forc_mean[max_i:max_f,::].collapsed('t',iris.analysis.MEAN)
TOf_min  = TO_forc_mean[min_i:min_f,::].collapsed('t',iris.analysis.MEAN)
TOr_max  = TO_rem_mean[max_i:max_f,::].collapsed('t',iris.analysis.MEAN)
TOr_min  = TO_rem_mean[min_i:min_f,::].collapsed('t',iris.analysis.MEAN)
TO_max  = TO_trop_mean[max_i:max_f,::].collapsed('t',iris.analysis.MEAN)
TO_min  = TO_trop_mean[min_i:min_f,::].collapsed('t',iris.analysis.MEAN)
TL_max  = TL_trop_mean[max_i:max_f,::].collapsed('t',iris.analysis.MEAN)
TL_min  = TL_trop_mean[min_i:min_f,::].collapsed('t',iris.analysis.MEAN)

plt.clf()
plt.ion()
plt.show()
plt.plot(TO_max.data,pressure,color='b')
plt.plot(TOf_max.data,pressure,color='r')
plt.plot(TOr_max.data,pressure,color='c')
# plt.plot(-TO_min.data,pressure,'--',color='b')
plt.plot(TL_max.data,pressure,color='g')
# plt.plot(-TL_min.data,pressure,'--',color='g')
# plt.legend((lat_range.astype(str)),loc=0,title='Latitude')
plt.title('Temp profile, Tropical Ocean, Max/Min forcing') 
plt.xlim(-0.2,1.1)
plt.ylabel('z [hPa]')
plt.xlabel('Temp')
plt.gca().invert_yaxis()


sys.exit('exiiiiiiiiit')

plt.figure(1)
qplt.pcmeshclf(Tmax,vmin=-1,vmax=1,cmap=mc.jetwhite())
plt.title('T response to max positive forcing')
plt.savefig('figures/comp_Tmax_sfc.png')

plt.figure(2)
qplt.pcmeshclf(Tmin,vmin=-1,vmax=1,cmap=mc.jetwhite())
plt.title('T response to max negative forcing')
plt.savefig('figures/comp_Tmin_sfc.png')

plt.figure(3)
qplt.pcmeshclf(T300max,vmin=-1,vmax=1,cmap=mc.jetwhite())
plt.title('T response to max positive forcing, 300hPa')
plt.savefig('figures/comp_Tmax_300.png')

plt.figure(4)
qplt.pcmeshclf(T300min,vmin=-1,vmax=1,cmap=mc.jetwhite())
plt.title('T response to max negative forcing, 300hPa')
plt.savefig('figures/comp_Tmin_300.png')

# plt.figure(3)
# qplt.pcmeshclf(T300max - T300min,vmin=-1,vmax=1,cmap=mc.jetwhite())
# plt.savefig('figures/reg_T_sfc_RH_sfc.png')









