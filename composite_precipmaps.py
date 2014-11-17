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

precip = iris.load_cube(ncfile_path + 'precip.4ysl.m48.nc')
rain_conv = iris.load_cube(ncfile_path + 'rain.conv.4ysl.m48.nc')
rain_lgscale = iris.load_cube(ncfile_path + 'rain.lgscale.4ysl.m48.nc')

# Define regions
precip.coord('latitude').guess_bounds()
precip.coord('longitude').guess_bounds()

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

precip_max            = precip[max_i:max_f,0,::].collapsed('time',iris.analysis.MEAN)
precip_min            = precip[min_i:min_f,0,::].collapsed('time',iris.analysis.MEAN)
rainconv_max            = rain_conv[max_i:max_f,0,::].collapsed('time',iris.analysis.MEAN)
rainconv_min            = rain_conv[min_i:min_f,0,::].collapsed('time',iris.analysis.MEAN)
rainlgsc_max            = rain_lgscale[max_i:max_f,0,::].collapsed('time',iris.analysis.MEAN)
rainlgsc_min            = rain_lgscale[min_i:min_f,0,::].collapsed('time',iris.analysis.MEAN)

rhmm = 1e-5
plt.figure(1)
qplt.pcmeshclf(precip_max,vmin=-rhmm,vmax=rhmm,cmap=mc.jetwhite())
plt.title('Precip response to max positive forcing. '+str(lag)+'-'+str(mons+lag)+' months')
plt.savefig('figures/comp_precip_max.png')

plt.figure(2)
qplt.pcmeshclf(precip_min,vmin=-rhmm,vmax=rhmm,cmap=mc.jetwhite_r())
plt.title('Precip response to min negative forcing. '+str(lag)+'-'+str(mons+lag)+' months')
plt.savefig('figures/comp_precip_min.png')

rhmm = 1e-5
plt.figure(3)
qplt.pcmeshclf(rainconv_max,vmin=-rhmm,vmax=rhmm,cmap=mc.jetwhite())
plt.title('Conv Rain response to max positive forcing. '+str(lag)+'-'+str(mons+lag)+' months')
plt.savefig('figures/comp_rainconv_max.png')

plt.figure(4)
qplt.pcmeshclf(rainconv_min,vmin=-rhmm,vmax=rhmm,cmap=mc.jetwhite_r())
plt.title('Conv Rain response to min negative forcing. '+str(lag)+'-'+str(mons+lag)+' months')
plt.savefig('figures/comp_rainconv_min.png')

rhmm = 2.5e-6
plt.figure(5)
qplt.pcmeshclf(rainlgsc_max,vmin=-rhmm,vmax=rhmm,cmap=mc.jetwhite())
plt.title('Largescale Rain response to max positive forcing. '+str(lag)+'-'+str(mons+lag)+' months')
plt.savefig('figures/comp_rainlgsc_max.png')

plt.figure(6)
qplt.pcmeshclf(rainlgsc_min,vmin=-rhmm,vmax=rhmm,cmap=mc.jetwhite_r())
plt.title('Largescale Rain response to min negative forcing. '+str(lag)+'-'+str(mons+lag)+' months')
plt.savefig('figures/comp_rainlgsc_min.png')



