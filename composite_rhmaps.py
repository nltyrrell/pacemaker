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


# Load in RHsfc, RHtropos (or RH300hPa?)
# calculate the regression between the timeseries
# RHsfc = alpha * RHtropo + epsilon

ncfile_path = '/home/nicholat/project/pacemaker/ncfiles/'
# rhum = iris.load_cube(ncfile_path + 'rhum.m48.sfc.4ysl.nc')
pres = iris.load_cube(ncfile_path + 'pres.sfc.4ysl.m48.nc')

# import the ACCESS data using iris
rhum_plv = iris.load_cube(ncfile_path + 'rhum.plv.4ysl.m48.nc')
rhum_plv.coord('p').standard_name = 'air_pressure'

# Define regions
rhum_plv.coord('latitude').guess_bounds()
rhum_plv.coord('longitude').guess_bounds()
pres.coord('latitude').guess_bounds()
pres.coord('longitude').guess_bounds()

# ---------- define sst/tland----------
lsmask = iris.load_cube(ncfile_path + 'lsmask.nc')[0,0,::]
landmask = ~(ma.make_mask(lsmask.data.copy()) + np.zeros(rhum_plv.shape)).astype(bool) # mask sea, show land
seamask = (ma.make_mask(lsmask.data.copy()) + np.zeros(rhum_plv.shape)).astype(bool) # mask land, show sea

# sys.exit('exiiiiiiitt')
RHocean = rhum_plv.copy()
RHland = rhum_plv.copy()
RHocean.data = ma.array(RHocean.data, mask=seamask)
RHland.data = ma.array(RHland.data, mask=landmask)
# --------------
tropics = iris.Constraint(latitude = lambda v: -30 <= v <= 30)
RHocean_trop = RHocean.extract(tropics)
grid_areas_trop = iris.analysis.cartography.area_weights(RHocean_trop)

mons = 6
lag = 0
max_i = 35 + lag
max_f = max_i + mons
min_i = 11 + lag
min_f = min_i + mons

# rhum_maxRH_mean  = rhum[max_i:max_f,0,::].collapsed('t',iris.analysis.MEAN)
# RHmax            = rhum_maxRH_mean
# rhum_minRH_mean  = rhum[min_i:min_f,0,::].collapsed('t',iris.analysis.MEAN)
# RHmin            = rhum_minRH_mean

rhum_300            = rhum_plv.extract(iris.Constraint(air_pressure=300))
RH300max            = rhum_300[max_i:max_f,::].collapsed('time',iris.analysis.MEAN)
RH300min            = rhum_300[min_i:min_f,::].collapsed('time',iris.analysis.MEAN)

rhum_700            = rhum_plv.extract(iris.Constraint(air_pressure=700))
RH700max            = rhum_700[max_i:max_f,::].collapsed('time',iris.analysis.MEAN)
RH700min            = rhum_700[min_i:min_f,::].collapsed('time',iris.analysis.MEAN)

rhum_1000            = rhum_plv.extract(iris.Constraint(air_pressure=1000))
RH1000max            = rhum_1000[max_i:max_f,::].collapsed('time',iris.analysis.MEAN)
RH1000min            = rhum_1000[min_i:min_f,::].collapsed('time',iris.analysis.MEAN)

Pmax  = pres[max_i:max_f,0,::].collapsed('time',iris.analysis.MEAN)
Pmin  = pres[min_i:min_f,0,::].collapsed('time',iris.analysis.MEAN)


# sys.exit('exit')
# plt.figure(1)
# qplt.pcmeshclf(RHmax,vmin=-1,vmax=1,cmap=mc.jetwhite())
# plt.title('RH response to max positive forcing')
# plt.savefig('figures/comp_RHmax_sfc.png')
# 
# plt.figure(2)
# qplt.pcmeshclf(RHmin,vmin=-1,vmax=1,cmap=mc.jetwhite_r())
# plt.title('RH response to min negative forcing')
# plt.savefig('figures/comp_RHmin_sfc.png')
# 
rhmm = 5
plt.figure(3)
qplt.pcmeshclf(RH300max,vmin=-rhmm,vmax=rhmm,cmap=mc.drywet2())
plt.title('RH response to max positive forcing, 300hPa. '+str(lag)+'-'+str(mons+lag)+' months')
plt.savefig('figures/comp_RHmax_300.pdf')

# plt.figure(4)
# qplt.pcmeshclf(RH300min,vmin=-rhmm,vmax=rhmm,cmap=mc.drywet_r())
# plt.title('RH response to min negative forcing, 300hPa. '+str(lag)+'-'+str(mons+lag)+' months')
# plt.savefig('figures/comp_RHmin_300.pdf')

plt.figure(5)
qplt.pcmeshclf(RH700max,vmin=-rhmm,vmax=rhmm,cmap=mc.drywet2())
plt.title('RH response to max positive forcing, 700hPa. '+str(lag)+'-'+str(mons+lag)+' months')
plt.savefig('figures/comp_RHmax_700.pdf')

# plt.figure(6)
# qplt.pcmeshclf(RH700min,vmin=-rhmm,vmax=rhmm,cmap=mc.drywet_r())
# plt.title('RH response to min negative forcing, 700hPa. '+str(lag)+'-'+str(mons+lag)+' months')
# plt.savefig('figures/comp_RHmin_700.pdf')

plt.figure(7)
qplt.pcmeshclf(RH1000max,vmin=-rhmm,vmax=rhmm,cmap=mc.drywet2())
plt.title('RH response to max positive forcing, 1000hPa. '+str(lag)+'-'+str(mons+lag)+' months')
plt.savefig('figures/comp_RHmax_1000.pdf')

# plt.figure(8)
# qplt.pcmeshclf(RH1000min,vmin=-rhmm,vmax=rhmm,cmap=mc.drywet_r())
# plt.title('RH response to min negative forcing, 1000hPa'+str(mons)+'-'+str(mons+lag)+' months')
# plt.savefig('figures/comp_RHmin_1000.pdf')

# plt.figure(3)
# qplt.pcmeshclf(RH300max - RH300min,vmin=-1,vmax=1,cmap=mc.jetwhite())
# plt.savefig('figures/reg_RH_sfc_RH_sfc.png')









