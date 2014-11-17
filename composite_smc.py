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

smc = iris.load_cube(ncfile_path + 'smc.sfc.4ysl.m48.nc')
smc_clim = iris.load_cube(ncfile_path + 'smc.sfc.4ysl.nc')

# Define regions
smc.coord('latitude').guess_bounds()
smc.coord('longitude').guess_bounds()

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

smc_max            = smc[max_i:max_f,0,::].collapsed('time',iris.analysis.MEAN)
smc_min            = smc[min_i:min_f,0,::].collapsed('time',iris.analysis.MEAN)
smc_mean            = smc_clim.collapsed('t',iris.analysis.MEAN)[0,::]


rhmm = 15
plt.figure(1)
qplt.pcmeshclf(smc_max,vmin=-rhmm,vmax=rhmm,cmap=mc.jetwhite())
plt.title('Soil moisture response to max positive forcing. '+str(lag)+'-'+str(mons+lag)+' months')
plt.savefig('figures/comp_smc_max.png')

plt.figure(2)
qplt.pcmeshclf(smc_min,vmin=-rhmm,vmax=rhmm,cmap=mc.jetwhite_r())
plt.title('Soil moisture response to min negative forcing. '+str(lag)+'-'+str(mons+lag)+' months')
plt.savefig('figures/comp_smc_min.png')

rhmm = 200
plt.figure(3)
qplt.pcmeshclf(smc_mean,vmin=-20,vmax=200,cmap=mc.drywet())
plt.title('Soil moisture climatology.')
plt.savefig('figures/smc_clim.png')

