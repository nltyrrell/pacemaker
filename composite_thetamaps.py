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
import ntiris as nti

# Load in Tsfc, Ttropos (or T300hPa?)
# calculate the regression between the timeseries
# Tsfc = alpha * Ttropo + epsilon

ncfile_path = '/home/nicholat/project/pacemaker/ncfiles/'
# theta = iris.load_cube(ncfile_path + 'rhum.m48.sfc.4ysl.nc')
# pres = iris.load_cube(ncfile_path + 'pres.sfc.4ysl.m48.nc')

# import the ACCESS data using iris
theta = iris.load_cube(ncfile_path + 'theta.plv.4ysl.m48.nc')
theta_anom = theta - theta.collapsed('time',iris.analysis.MEAN)
theta_plv = nti.remove_seascyc(theta_anom)
thetav = iris.load_cube(ncfile_path + 'thetav.plv.4ysl.m48.nc')
thetav_anom = thetav - thetav.collapsed('time',iris.analysis.MEAN)
thetav_plv = nti.remove_seascyc(thetav_anom)
temp = iris.load_cube(ncfile_path + 'temp.plv.4ysl.m48.nc')
temp_anom = temp - temp.collapsed('time',iris.analysis.MEAN)
temp_plv = nti.remove_seascyc(temp_anom)
try:
    theta_plv.coord('t').standard_name='time'
except:
    pass
else:
    print "t coord changed to time"

# Define regions
theta_plv.coord('latitude').guess_bounds()
theta_plv.coord('longitude').guess_bounds()
thetav_plv.coord('latitude').guess_bounds()
thetav_plv.coord('longitude').guess_bounds()
temp_plv.coord('latitude').guess_bounds()
temp_plv.coord('longitude').guess_bounds()
# pres.coord('latitude').guess_bounds()
# pres.coord('longitude').guess_bounds()

# ---------- define sst/tland----------
lsmask = iris.load_cube(ncfile_path + 'lsmask.nc')[0,0,::]
landmask = ~(ma.make_mask(lsmask.data.copy()) + np.zeros(theta_plv.shape)).astype(bool) # mask sea, show land
seamask = (ma.make_mask(lsmask.data.copy()) + np.zeros(theta_plv.shape)).astype(bool) # mask land, show sea

# sys.exit('exiiiiiiitt')
Tocean = theta_plv.copy()
Tland = theta_plv.copy()
Tocean.data = ma.array(Tocean.data, mask=seamask)
Tland.data = ma.array(Tland.data, mask=landmask)
# --------------
tropics = iris.Constraint(latitude = lambda v: -30 <= v <= 30)
Tocean_trop = Tocean.extract(tropics)
grid_areas_trop = iris.analysis.cartography.area_weights(Tocean_trop)

mons = 3
lag = 3
max_i = 35 + lag
max_f = max_i + mons
min_i = 11 + lag
min_f = min_i + mons

# theta_maxT_mean  = theta[max_i:max_f,0,::].collapsed('t',iris.analysis.MEAN)
# Tmax            = theta_maxT_mean
# theta_minT_mean  = theta[min_i:min_f,0,::].collapsed('t',iris.analysis.MEAN)
# Tmin            = theta_minT_mean

theta_300            = theta_plv.extract(iris.Constraint(air_pressure=300))
Th300max  = theta_300[max_i:max_f,::].collapsed('time',iris.analysis.MEAN)
Th300min  = theta_300[min_i:min_f,::].collapsed('time',iris.analysis.MEAN)

theta_700            = theta_plv.extract(iris.Constraint(air_pressure=700))
Th700max  = theta_700[max_i:max_f,::].collapsed('time',iris.analysis.MEAN)
Th700min  = theta_700[min_i:min_f,::].collapsed('time',iris.analysis.MEAN)

theta_300            = thetav_plv.extract(iris.Constraint(air_pressure=300))
Thv300max  = theta_300[max_i:max_f,::].collapsed('time',iris.analysis.MEAN)
Thv300min  = theta_300[min_i:min_f,::].collapsed('time',iris.analysis.MEAN)

theta_700            = thetav_plv.extract(iris.Constraint(air_pressure=700))
Thv700max  = theta_700[max_i:max_f,::].collapsed('time',iris.analysis.MEAN)
Thv700min  = theta_700[min_i:min_f,::].collapsed('time',iris.analysis.MEAN)

temp_300            = theta_plv.extract(iris.Constraint(air_pressure=300))
temp_300_maxT_mean  = theta_300[max_i:max_f,::].collapsed('time',iris.analysis.MEAN)
T300max             = temp_300_maxT_mean
temp_300_minT_mean  = theta_300[min_i:min_f,::].collapsed('time',iris.analysis.MEAN)
T300min             = temp_300_minT_mean

# Pmax  = pres[max_i:max_f,0,::].collapsed('time',iris.analysis.MEAN)
# Pmin  = pres[min_i:min_f,0,::].collapsed('time',iris.analysis.MEAN)


# sys.exit('exit')
# plt.figure(1)
# qplt.pcmeshclf(Tmax,vmin=-1,vmax=1,cmap=mc.jetwhite())
# plt.title('T response to max positive forcing')
# plt.savefig('figures/comp_Tmax_sfc.png')
# 
# plt.figure(2)
# qplt.pcmeshclf(Tmin,vmin=-1,vmax=1,cmap=mc.jetwhite_r())
# plt.title('T response to min negative forcing')
# plt.savefig('figures/comp_Tmin_sfc.png')
# 
thmm = 2
plt.close('all')

plt.figure(1)
qplt.pcmeshclf(Th300max,vmin=-thmm,vmax=thmm,cmap=mc.jetwhite())
plt.title('Theta response to max positive forcing, 300hPa')
# plt.savefig('figures/comp_theta_max_300.png')

plt.figure(2)
qplt.pcmeshclf(Th300min,vmin=-thmm,vmax=thmm,cmap=mc.jetwhite_r())
plt.title('Theta response to min negative forcing, 300hPa')
# plt.savefig('figures/comp_theta_min_300.png')

plt.figure(3)
qplt.pcmeshclf(Th700max,vmin=-thmm,vmax=thmm,cmap=mc.jetwhite())
plt.title('Theta response to max positive forcing, 700hPa')
# plt.savefig('figures/comp_theta_max_700.png')

plt.figure(4)
qplt.pcmeshclf(Th700min,vmin=-thmm,vmax=thmm,cmap=mc.jetwhite_r())
plt.title('Theta response to min negative forcing, 700hPa')
# plt.savefig('figures/comp_theta_min_700.png')

plt.figure(5)
qplt.pcmeshclf(Thv300max,vmin=-thmm,vmax=thmm,cmap=mc.jetwhite())
plt.title('Theta v response to max positive forcing, 300hPa')
# plt.savefig('figures/comp_theta_max_300.png')

plt.figure(6)
qplt.pcmeshclf(Thv300min,vmin=-thmm,vmax=thmm,cmap=mc.jetwhite_r())
plt.title('Theta v response to min negative forcing, 300hPa')
# plt.savefig('figures/comp_theta_min_300.png')

plt.figure(7)
qplt.pcmeshclf(Thv700max,vmin=-thmm,vmax=thmm,cmap=mc.jetwhite())
plt.title('Theta v response to max positive forcing, 700hPa')
# plt.savefig('figures/comp_theta_max_700.png')

plt.figure(8)
qplt.pcmeshclf(Thv700min,vmin=-thmm,vmax=thmm,cmap=mc.jetwhite_r())
plt.title('Theta v response to min negative forcing, 700hPa')
# plt.savefig('figures/comp_theta_min_700.png')

# plt.figure(3)
# qplt.pcmeshclf(T300max - T300min,vmin=-1,vmax=1,cmap=mc.jetwhite())
# plt.savefig('figures/reg_T_sfc_RH_sfc.png')









