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


# Load in Gsfc, Gtropos (or G300hPa?)
# calculate the regression between the timeseries
# Gsfc = alpha * Gtropo + epsilon

ncfile_path = '/home/nicholat/project/pacemaker/ncfiles/'
# gpht = iris.load_cube(ncfile_path + 'rhum.m48.sfc.4ysl.nc')
gpht = iris.load_cube(ncfile_path + 'gpht.plv.4ysl.m48.nc')

# import the ACCESS data using iris
gpht_plv = iris.load_cube(ncfile_path + 'gpht.plv.4ysl.m48.nc')
gpht_plv.coord('p').standard_name = 'air_pressure'

# Define regions
gpht_plv.coord('latitude').guess_bounds()
gpht_plv.coord('longitude').guess_bounds()
# pres.coord('latitude').guess_bounds()
# pres.coord('longitude').guess_bounds()

# ---------- define sst/tland----------
lsmask = iris.load_cube(ncfile_path + 'lsmask.nc')[0,0,::]
landmask = ~(ma.make_mask(lsmask.data.copy()) + np.zeros(gpht_plv.shape)).astype(bool) # mask sea, show land
seamask = (ma.make_mask(lsmask.data.copy()) + np.zeros(gpht_plv.shape)).astype(bool) # mask land, show sea

# sys.exit('exiiiiiiitt')
Gocean = gpht_plv.copy()
Gland = gpht_plv.copy()
Gocean.data = ma.array(Gocean.data, mask=seamask)
Gland.data = ma.array(Gland.data, mask=landmask)
# --------------
tropics = iris.Constraint(latitude = lambda v: -30 <= v <= 30)
Gocean_trop = Gocean.extract(tropics)
grid_areas_trop = iris.analysis.cartography.area_weights(Gocean_trop)

mons = 3
lag = 3
max_i = 35 + lag
max_f = max_i + mons
min_i = 11 + lag
min_f = min_i + mons

# gpht_maxG_mean  = gpht[max_i:max_f,0,::].collapsed('t',iris.analysis.MEAN)
# Gmax            = gpht_maxG_mean
# gpht_minG_mean  = gpht[min_i:min_f,0,::].collapsed('t',iris.analysis.MEAN)
# Gmin            = gpht_minG_mean

gpht_500            = gpht_plv.extract(iris.Constraint(air_pressure=500))
gpht_500_maxG_mean  = gpht_500[max_i:max_f,::].collapsed('time',iris.analysis.MEAN)
G500max             = gpht_500_maxG_mean
gpht_500_minG_mean  = gpht_500[min_i:min_f,::].collapsed('time',iris.analysis.MEAN)
G500min             = gpht_500_minG_mean

gpht_300            = gpht_plv.extract(iris.Constraint(air_pressure=300))
gpht_300_maxG_mean  = gpht_300[max_i:max_f,::].collapsed('time',iris.analysis.MEAN)
G300max             = gpht_300_maxG_mean
gpht_300_minG_mean  = gpht_300[min_i:min_f,::].collapsed('time',iris.analysis.MEAN)
G300min             = gpht_300_minG_mean

gpht_700            = gpht_plv.extract(iris.Constraint(air_pressure=700))
gpht_700_maxG_mean  = gpht_700[max_i:max_f,::].collapsed('time',iris.analysis.MEAN)
G700max             = gpht_700_maxG_mean
gpht_700_minG_mean  = gpht_700[min_i:min_f,::].collapsed('time',iris.analysis.MEAN)
G700min             = gpht_700_minG_mean

# Pmax  = pres[max_i:max_f,0,::].collapsed('time',iris.analysis.MEAN)
# Pmin  = pres[min_i:min_f,0,::].collapsed('time',iris.analysis.MEAN)


# sys.exit('exit')
# plt.figure(1)
# qplt.pcmeshclf(Gmax,vmin=-1,vmax=1,cmap=mc.jetwhite())
# plt.title('G response to max positive forcing')
# plt.savefig('figures/comp_Gmax_sfc.png')
# 
# plt.figure(2)
# qplt.pcmeshclf(Gmin,vmin=-1,vmax=1,cmap=mc.jetwhite_r())
# plt.title('G response to min negative forcing')
# plt.savefig('figures/comp_Gmin_sfc.png')
# 
gpmm = 50
plt.figure(3)
qplt.pcmeshclf(G500max,vmin=-gpmm,vmax=gpmm,cmap=mc.jetwhite())
plt.title('GPHT response to max positive forcing, 500hPa')
plt.savefig('figures/comp_GPHTmax_500.pdf')

plt.figure(4)
qplt.pcmeshclf(G300max,vmin=-gpmm,vmax=gpmm,cmap=mc.jetwhite_r())
plt.title('GPHT response to max negative forcing, 300hPa')
plt.savefig('figures/comp_GPHTmax_300.png')

plt.figure(5)
qplt.pcmeshclf(G700max,vmin=-gpmm,vmax=gpmm,cmap=mc.jetwhite())
plt.title('GPHT response to max positive forcing, 700hPa')
plt.savefig('figures/comp_GPHTmax_700.png')

sys.exit()
plt.figure(6)
qplt.pcmeshclf(G700min,vmin=-gpmm,vmax=gpmm,cmap=mc.jetwhite_r())
plt.title('GPHT response to min negative forcing, 700hPa')
plt.savefig('figures/comp_GPHTmin_700.png')

# plt.figure(3)
# qplt.pcmeshclf(G300max - G300min,vmin=-1,vmax=1,cmap=mc.jetwhite())
# plt.savefig('figures/reg_G_sfc_RH_sfc.png')









