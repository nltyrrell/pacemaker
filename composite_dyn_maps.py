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
cld = iris.load_cube(ncfile_path + 'cld.thlev.4ysl.m48.nc')
slp = iris.load_cube(ncfile_path + 'pres.sfc.4ysl.m48.nc')
cld.coord('Hybrid height').standard_name = 'atmosphere_hybrid_height_coordinate'

# Define regions
cld.coord('latitude').guess_bounds()
cld.coord('longitude').guess_bounds()
slp.coord('latitude').guess_bounds()
slp.coord('longitude').guess_bounds()

high = iris.Constraint(atmosphere_hybrid_height_coordinate = lambda h: 5000 <= h <= 15000)
low = iris.Constraint(atmosphere_hybrid_height_coordinate = lambda h: 0 <= h <= 5000)
full = iris.Constraint(atmosphere_hybrid_height_coordinate = lambda h: 0 <= h <= 10000)

cld_high = cld.extract(high).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN).data
cld_low = cld.extract(low).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN).data
cld_full = cld.extract(full).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN).data

mons = 3
lag = 3
max_i = 35 + lag
max_f = max_i + mons
min_i = 11 + lag
min_f = min_i + mons

slp_max = slp[max_i:max_f,0,::].collapsed('time',iris.analysis.MEAN)
slp_min = slp[min_i:min_f,0,::].collapsed('time',iris.analysis.MEAN)
cld_high_max = cld_high[max_i:max_f,::].mean(axis=0)
cld_high_min = cld_high[min_i:min_f,::].mean(axis=0)
cld_low_max = cld_low[max_i:max_f,::].mean(axis=0)
cld_low_min = cld_low[min_i:min_f,::].mean(axis=0)
cld_full_max = cld_full[max_i:max_f,::].mean(axis=0)
cld_full_min = cld_full[min_i:min_f,::].mean(axis=0)

cld_high_max_cube = slp_max.copy()
cld_high_min_cube = slp_max.copy()
cld_low_max_cube = slp_max.copy()
cld_low_min_cube = slp_max.copy()
cld_full_max_cube = slp_max.copy()
cld_full_min_cube = slp_max.copy()

cld_high_max_cube.data = cld_high_max
cld_high_min_cube.data = cld_high_min
cld_low_max_cube.data = cld_low_max
cld_low_min_cube.data = cld_low_min
cld_full_max_cube.data = cld_full_max
cld_full_min_cube.data = cld_full_min

plt.clf()
plt.ion()

plt.figure(1)
qplt.pcmeshclf(slp_max,vmin=-270,vmax=270,cmap=mc.jetwhite())
plt.title('SLP response to max positive forcing')
# plt.savefig('figures/comp_slpmax.png')

plt.figure(2)
qplt.pcmeshclf(slp_min,vmin=-270,vmax=270,cmap=mc.jetwhite_r())
plt.title('SLP response to min positive forcing')
# plt.savefig('figures/comp_slpmin.png')

plt.figure(3)
qplt.pcmeshclf(cld_low_max_cube,vmin=-0.025,vmax=0.025,cmap=mc.jetwhite())
plt.title('Low cloud response to max positive forcing')
# plt.savefig('figures/comp_cld_lowmax.png')

plt.figure(4)
qplt.pcmeshclf(cld_low_min_cube,vmin=-0.025,vmax=0.025,cmap=mc.jetwhite_r())
plt.title('Low cloud response to min positive forcing')
# plt.savefig('figures/comp_cld_lowmin.png')

plt.figure(5)
qplt.pcmeshclf(cld_high_max_cube,vmin=-0.025,vmax=0.025,cmap=mc.jetwhite())
plt.title('High cloud response to max positive forcing')
# plt.savefig('figures/comp_cld_highmax.png')

plt.figure(6)
qplt.pcmeshclf(cld_high_min_cube,vmin=-0.025,vmax=0.025,cmap=mc.jetwhite_r())
plt.title('High cloud response to min positive forcing')
# plt.savefig('figures/comp_cld_lowmin.png')









