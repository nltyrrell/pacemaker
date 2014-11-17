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

lsclsk = iris.load_cube(ncfile_path + 'lwflux.clsky.sfc.4ysl.m48.nc')
dlwr = iris.load_cube(ncfile_path + 'dlwr.sfc.4ysl.m48.nc')
dswr = iris.load_cube(ncfile_path + 'dswr.sfc.4ysl.m48.nc')

# Define regions
lsclsk.coord('latitude').guess_bounds()
lsclsk.coord('longitude').guess_bounds()
dlwr.coord('latitude').guess_bounds()
dlwr.coord('longitude').guess_bounds()
dswr.coord('latitude').guess_bounds()
dswr.coord('longitude').guess_bounds()

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

lwclsk_max            = lsclsk[max_i:max_f,0,::].collapsed('time',iris.analysis.MEAN)
lwclsk_min            = lsclsk[min_i:min_f,0,::].collapsed('time',iris.analysis.MEAN)
dlwr_max            = dlwr[max_i:max_f,0,::].collapsed('time',iris.analysis.MEAN)
dlwr_min            = dlwr[min_i:min_f,0,::].collapsed('time',iris.analysis.MEAN)
dswr_max            = dswr[max_i:max_f,0,::].collapsed('time',iris.analysis.MEAN)
dswr_min            = dswr[min_i:min_f,0,::].collapsed('time',iris.analysis.MEAN)


rhmm = 10
plt.figure(1)
qplt.pcmeshclf(lwclsk_max,vmin=-rhmm,vmax=rhmm,cmap=mc.jetwhite())
plt.title('LW clearsky response to max positive forcing. '+str(lag)+'-'+str(mons+lag)+' months')
plt.savefig('figures/comp_lwclsky_max.png')

plt.figure(2)
qplt.pcmeshclf(lwclsk_min,vmin=-rhmm,vmax=rhmm,cmap=mc.jetwhite_r())
plt.title('LW clearsky response to min negative forcing. '+str(lag)+'-'+str(mons+lag)+' months')
plt.savefig('figures/comp_lwclsky_min.png')

plt.figure(3)
qplt.pcmeshclf(dlwr_max,vmin=-rhmm,vmax=rhmm,cmap=mc.jetwhite())
plt.title('Down LW response to max positive forcing. '+str(lag)+'-'+str(mons+lag)+' months')
plt.savefig('figures/comp_dlwr_max.png')

plt.figure(4)
qplt.pcmeshclf(dlwr_min,vmin=-rhmm,vmax=rhmm,cmap=mc.jetwhite_r())
plt.title('Down LW response to min negative forcing. '+str(lag)+'-'+str(mons+lag)+' months')
plt.savefig('figures/comp_dlwr_min.png')

plt.figure(5)
qplt.pcmeshclf(dswr_max,vmin=-rhmm,vmax=rhmm,cmap=mc.jetwhite())
plt.title('Down SW response to max positive forcing. '+str(lag)+'-'+str(mons+lag)+' months')
plt.savefig('figures/comp_dswr_max.png')

plt.figure(6)
qplt.pcmeshclf(dswr_min,vmin=-rhmm,vmax=rhmm,cmap=mc.jetwhite_r())
plt.title('Down SW response to min negative forcing. '+str(lag)+'-'+str(mons+lag)+' months')
plt.savefig('figures/comp_dswr_min.png')


