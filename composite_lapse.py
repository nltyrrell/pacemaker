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


ncfile_path = '/home/nicholat/project/pacemaker/ncfiles/'
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

trop_forcing = iris.Constraint(latitude = lambda v: -20 <= v <= 20,longitude = lambda l: 140 <= l <= 300)
Tocean_forc = Tocean.extract(trop_forcing)
grid_areas_forc = iris.analysis.cartography.area_weights(Tocean_forc)
TO_forc_mean = Tocean_forc.collapsed(['latitude','longitude'],
                iris.analysis.MEAN,weights=grid_areas_forc)

trop_remote = iris.Constraint(latitude = lambda v: -30 <= v <= 30,
                longitude = lambda l: (300 <= l <= 360) or (0 <= l <= 130) )
Tocean_rem = Tocean.extract(trop_remote)
grid_areas_rem = iris.analysis.cartography.area_weights(Tocean_rem)
TO_rem_mean = Tocean_rem.collapsed(['latitude','longitude'],
                iris.analysis.MEAN,weights=grid_areas_rem)

Tland_trop = Tland.extract(tropics)
grid_areas_trop = iris.analysis.cartography.area_weights(Tland_trop)
TL_trop_mean = Tland_trop.collapsed(['latitude','longitude'],
                iris.analysis.MEAN,weights=grid_areas_trop)
mons = 6
lag = 0
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
# plt.plot(TO_max.data,pressure,color='b')
# plt.plot(-TO_min.data,pressure,'--',color='b',label="_nolegend_")
plt.plot(TOf_max.data,pressure,color='r',linewidth=1.8,label='Tocean forcing')
plt.plot(-TOf_min.data,pressure,'--',color='r',linewidth=1.8,label="_nolegend_")
plt.plot(TOr_max.data,pressure,color='c',linewidth=1.8,label='Tocean remote')
plt.plot(-TOr_min.data,pressure,'--',color='c',label="_nolegend_",linewidth=1.8)
plt.plot(TL_max.data,pressure,color='g',linewidth=1.8,label='Tland')
plt.plot(-TL_min.data,pressure,'--',color='g',label="_nolegend_",linewidth=1.8)
ext1 = plt.plot(0,0,'-',linewidth=1,color='k',label='Max')
ext2 = plt.plot(0,0,'--',linewidth=1,color='k',label='Min')
plt.legend(loc=0,title='Latitude')
plt.title('Tropics, Max/Min forcing, '+str(lag)+' to '+str(lag+mons)+' months from peak forcing') 
plt.xlim(-0.2,1.5)
plt.ylabel('z [hPa]')
plt.xlabel('Temp')
plt.gca().invert_yaxis()
plt.savefig('./figures/comp_trop_zonal_profiles_m'+str(lag)+str(lag+mons)+'.png')

