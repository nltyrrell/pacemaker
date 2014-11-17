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

# import the ACCESS data using iris
cld = iris.load_cube(ncfile_path + 'cld.thlev.4ysl.m48.nc')
cld.coord('Hybrid height').standard_name = 'atmosphere_hybrid_height_coordinate'

# Define regions
cld.coord('latitude').guess_bounds()
cld.coord('longitude').guess_bounds()

mons = 6
lag = 0
max_i = 35 + lag
max_f = max_i + mons
min_i = 11 + lag
min_f = min_i + mons

high = iris.Constraint(atmosphere_hybrid_height_coordinate = lambda h: 5000 <= h <= 15000)
low = iris.Constraint(atmosphere_hybrid_height_coordinate = lambda h: 0 <= h <= 5000)
full = iris.Constraint(atmosphere_hybrid_height_coordinate = lambda h: 0 <= h <= 10000)

cld_max            = cld[max_i:max_f,::].collapsed('time',iris.analysis.MEAN)
cld_min            = cld[min_i:min_f,::].collapsed('time',iris.analysis.MEAN)

cld_high_max = cld_max.extract(high).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)
cld_high_min = cld_min.extract(high).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)

cld_low_max = cld_max.extract(low).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)
cld_low_min = cld_min.extract(low).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)

cld_full_max = cld_max.extract(full).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)
cld_full_min = cld_min.extract(full).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)

rhmm = 0.02
plt.figure(3)
qplt.pcmeshclf(cld_high_max,vmin=-rhmm,vmax=rhmm,cmap=mc.jetwhite())
plt.title('High Cloud response to max positive forcing. '+str(lag)+'-'+str(mons+lag)+' months')
plt.savefig('figures/comp_cldhigh_max.png')

plt.figure(4)
qplt.pcmeshclf(cld_high_min,vmin=-rhmm,vmax=rhmm,cmap=mc.jetwhite_r())
plt.title('High Cloud response to min negative forcing. '+str(lag)+'-'+str(mons+lag)+' months')
plt.savefig('figures/comp_cldhigh_min.png')

plt.figure(5)
qplt.pcmeshclf(cld_low_max,vmin=-rhmm,vmax=rhmm,cmap=mc.jetwhite())
plt.title('Low Cloud response to max positive forcing, 700hPa. '+str(lag)+'-'+str(mons+lag)+' months')
plt.savefig('figures/comp_cldlow_max.png')

plt.figure(6)
qplt.pcmeshclf(cld_low_min,vmin=-rhmm,vmax=rhmm,cmap=mc.jetwhite_r())
plt.title('Low Cloud response to min negative forcing, 700hPa. '+str(lag)+'-'+str(mons+lag)+' months')
plt.savefig('figures/comp_cldlow_min.png')

plt.figure(7)
qplt.pcmeshclf(cld_full_max,vmin=-rhmm,vmax=rhmm,cmap=mc.jetwhite())
plt.title('All Cloud response to max positive forcing, 1000hPa. '+str(lag)+'-'+str(mons+lag)+' months')
plt.savefig('figures/comp_cldfull_max.png')

plt.figure(8)
qplt.pcmeshclf(cld_full_min,vmin=-rhmm,vmax=rhmm,cmap=mc.jetwhite_r())
plt.title('All Cloud response to min negative forcing, 1000hPa'+str(mons)+'-'+str(mons+lag)+' months')
plt.savefig('figures/comp_cldfull_min.png')

# plt.figure(3)
# qplt.pcmeshclf(RH300max - RH300min,vmin=-1,vmax=1,cmap=mc.jetwhite())
# plt.savefig('figures/reg_RH_sfc_RH_sfc.png')









