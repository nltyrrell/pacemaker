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
temp = iris.load_cube(ncfile_path + 'temp.sfc.4ysl.ym.nc')
temp_mean = temp[:,0,::].collapsed('time',iris.analysis.MEAN)

# import the ACCESS data using iris
temp_plv = iris.load_cube(ncfile_path + 'temp.m48.plv.4ysl.nc')
temp_plv.coord('p').standard_name = 'air_pressure'

temp_plv_mean = temp_plv[:,:,::].collapsed('t',iris.analysis.MEAN)
temp_plv_anom = temp_plv-temp_plv_mean

# Define regions
temp_plv_anom.coord('latitude').guess_bounds()
temp_plv_anom.coord('longitude').guess_bounds()

# ---------- define sst/tland----------
lsmask = iris.load_cube(ncfile_path + 'lsmask.nc')[0,0,::]
landmask = ~(ma.make_mask(lsmask.data.copy()) + np.zeros(temp_plv.shape)).astype(bool) # mask sea, show land
seamask = (ma.make_mask(lsmask.data.copy()) + np.zeros(temp_plv.shape)).astype(bool) # mask land, show sea

Tocean = temp_plv_anom.copy()
Tland = temp_plv_anom.copy()
Tocean.data = ma.array(Tocean.data, mask=seamask)
Tland.data = ma.array(Tland.data, mask=landmask)
# --------------

tropics = iris.Constraint(latitude = lambda v: -30 <= v <= 30)
Tocean_trop = Tocean.extract(tropics)
grid_areas_trop = iris.analysis.cartography.area_weights(Tocean_trop)

Tocean_trop_mean = Tocean_trop.collapsed(['latitude', 'longitude'],
                           iris.analysis.MEAN,
                           weights=grid_areas_trop)


SthAm = iris.Constraint(longitude=lambda l: (290 <= l <= 315), latitude = lambda l: (-23 <= l <= 0))
Tplv_SthAm = temp_plv_anom.extract(SthAm)
grid_areas_SthAm = iris.analysis.cartography.area_weights(Tplv_SthAm)

Tplv_SthAm_armean = Tplv_SthAm.collapsed(['latitude', 'longitude'],
                           iris.analysis.MEAN,
                           weights=grid_areas_SthAm)

NthAfr = iris.Constraint(longitude=lambda l: (0 <= l <= 25), latitude = lambda l: (10 <= l <= 30))
Tplv_NthAfr = temp_plv_anom.extract(NthAfr)
grid_areas_NthAfr = iris.analysis.cartography.area_weights(Tplv_NthAfr)

Tplv_NthAfr_armean = Tplv_NthAfr.collapsed(['latitude', 'longitude'],
                           iris.analysis.MEAN,
                           weights=grid_areas_NthAfr)

CntAfr = iris.Constraint(longitude=lambda l: (12 <= l <= 40), latitude = lambda l: (-15 <= l <= 5))
Tplv_CntAfr = temp_plv_anom.extract(CntAfr)
grid_areas_CntAfr = iris.analysis.cartography.area_weights(Tplv_CntAfr)

Tplv_CntAfr_armean = Tplv_CntAfr.collapsed(['latitude', 'longitude'],
                           iris.analysis.MEAN,
                           weights=grid_areas_CntAfr)

Aus = iris.Constraint(longitude=lambda l: (120 <= l <= 140), latitude = lambda l: (-30 <= l <= -17))
Tplv_Aus = temp_plv_anom.extract(Aus)
grid_areas_Aus = iris.analysis.cartography.area_weights(Tplv_Aus)

Tplv_Aus_armean = Tplv_Aus.collapsed(['latitude', 'longitude'],
                           iris.analysis.MEAN,
                           weights=grid_areas_Aus)

hold = {}
# temp[height (100), evfr (4), ssts (11), l/o col (0-Oc,1-Ln)
# col which is used for regresion, '0' for ocean, '1' for land
lo_col = 0

NthAfrhold = np.zeros((2,temp_plv.coord('air_pressure').shape[0]))
CntAfrhold = np.zeros((2,temp_plv.coord('air_pressure').shape[0]))
SthAmhold = np.zeros((2,temp_plv.coord('air_pressure').shape[0]))
Aushold = np.zeros((2,temp_plv.coord('air_pressure').shape[0]))
TrTochold = np.zeros(temp_plv.coord('air_pressure').shape[0])


mons = 6
lag = 0
max_i = 35 + lag
max_f = max_i + mons
min_i = 11 + lag
min_f = min_i + mons
TOMax = Tocean_trop_mean[max_i:max_f,::].collapsed('t',iris.analysis.MEAN)
SthAmMax = Tplv_SthAm_armean[max_i:max_f,::].collapsed('t',iris.analysis.MEAN)
NthAfrMax = Tplv_NthAfr_armean[max_i:max_f,::].collapsed('t',iris.analysis.MEAN)
CntAfrMax = Tplv_CntAfr_armean[max_i:max_f,::].collapsed('t',iris.analysis.MEAN)
AusMax = Tplv_Aus_armean[max_i:max_f,::].collapsed('t',iris.analysis.MEAN)

pressure = temp_plv.coord('air_pressure').points
plt.clf()
plt.plot(NthAfrMax.data,pressure,'g')#linestyle='--',color=wetdry(k,nloops))
plt.plot(CntAfrMax.data,pressure,'b')#linestyle='--',color=wetdry(k,nloops))
plt.plot(SthAmMax.data,pressure,'r')#linestyle='--',color=wetdry(k,nloops))
plt.plot(AusMax.data,pressure,'k')#linestyle='--',color=wetdry(k,nloops))
plt.plot(TOMax.data,pressure,'c')#linestyle='--',color=wetdry(k,nloops))

plt.legend(('NthAfr','Central Africa','Sth America','Australia','Tropical Ocean'),loc=3, )
plt.title('Temp profile, time mean for 0 to 6 months from max forcing') 
plt.xlim(-0.5,0.7)
plt.ylabel('z [hPa]')
plt.xlabel('Regresion coeff.')
plt.gca().invert_yaxis()
plt.savefig('./figures/cont_maxtemp_profile.png')




