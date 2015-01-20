import numpy as np
import numpy.ma as ma
import iris as iris
import iris.plot as iplt
import iris.quickplot as qplt
from windspharm.iris import VectorWind
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import mycmaps as mc
import scipy.stats as stats
import sys
import troposave as ta
import mapping as mp
import ntiris as nt

ncfile_path = '/home/nicholat/project/pacemaker/ncfiles/'
w = iris.load_cube(ncfile_path + 'w.thlev.4ysl.nc')
tsfc = iris.load_cube(ncfile_path + 'temp.sfc.4ysl.nc')
w.coords()[0].standard_name = 'time'
tsfc.coords()[0].standard_name = 'time'

thlev = w.coord('atmosphere_hybrid_height_coordinate').points

tsfc = nt.remove_seascyc(tsfc)
w = nt.remove_seascyc(w)

# Define regions
w.coord('latitude').guess_bounds()
w.coord('longitude').guess_bounds()
tsfc.coord('latitude').guess_bounds()
tsfc.coord('longitude').guess_bounds()

# ---------- define sst/tland----------
lsmask = iris.load_cube(ncfile_path + 'lsmask.nc')[0,0,::]
lsmask.coord('latitude').guess_bounds()
lsmask.coord('longitude').guess_bounds()
landmask = ~(ma.make_mask(lsmask.data.copy()) + np.zeros(tsfc.shape)).astype(bool) # mask sea, show land
seamask = (ma.make_mask(lsmask.data.copy()) + np.zeros(tsfc.shape)).astype(bool) # mask land, show sea

Tocean = tsfc.copy()
Tland = tsfc.copy()
Tocean.data = ma.array(Tocean.data, mask=seamask)
Tland.data = ma.array(Tland.data, mask=landmask)
# landmask = ~(ma.make_mask(lsmask.data.copy()) + np.zeros(temp_plv.shape)).astype(bool) # mask sea, show land
# seamask = (ma.make_mask(lsmask.data.copy()) + np.zeros(temp_plv.shape)).astype(bool) # mask land, show sea

# define upper, middle, lower regions of the troposphere for u, v
high = iris.Constraint(atmosphere_hybrid_height_coordinate = lambda h: 5000 <= h <= 15000)
low = iris.Constraint(atmosphere_hybrid_height_coordinate = lambda h: 0 <= h <= 5000)
w_up = w.extract(high).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)
w_low = w.extract(low).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)
 
# Define Australian region and generate timeseries of mean Aus sfc temps
temp_sfc_anom = Tland - Tland.collapsed('time',iris.analysis.MEAN)
Aus = iris.Constraint(longitude=lambda l: (109 <= l <= 160), latitude = lambda l: (-42 <= l <= -11))
Tsfc_Aus = temp_sfc_anom.extract(Aus)
grid_areas_Aus = iris.analysis.cartography.area_weights(Tsfc_Aus)

Tsfc_Aus_mean = Tsfc_Aus.collapsed(['latitude', 'longitude'],
                           iris.analysis.MEAN,
                           weights=grid_areas_Aus)


regw_low, corw_low,  pvalw_low = nt.linregts(w_low,Tsfc_Aus_mean[:,0],'w_lowtrop','tsfc_Aus')
# regw_mid, corw_mid,  pvalw_mid = nt.linregts(w_midtrop,Tsfc_Aus_mean[:,0],'w_midtrop','tsfc_Aus')
regw_up,  corw_up,  pvalw_up  = nt.linregts(w_up,Tsfc_Aus_mean[:,0],'w_uptrop','tsfc_Aus')

plt.figure(1)
qplt.pcmeshclf(regw_low,vmin=-100,vmax=100,cmap=mc.jetwhite())
plt.title('Regression between Tsfc Aus, w 0-5000m')
plt.savefig('./figures/aus_wlow_reg.pdf')
plt.figure(2)
qplt.pcmeshclf(corw_low,vmin=-0.7,vmax=0.7,cmap=mc.jetwhite())
plt.title('Correlation between Tsfc Aus, w 0-5000m')
plt.savefig('./figures/aus_wlow_cor.pdf')
plt.figure(3)
qplt.pcmeshclf(regw_up,vmin=-100,vmax=100,cmap=mc.jetwhite())
plt.title('Regression between Tsfc Aus, w 5000-15000m')
plt.savefig('./figures/aus_wup_reg.pdf')
plt.figure(4)
qplt.pcmeshclf(corw_up,vmin=-0.7,vmax=0.7,cmap=mc.jetwhite())
plt.title('Correlation between Tsfc Aus, w 5000-15000m')
plt.savefig('./figures/aus_wup_cor.pdf')






