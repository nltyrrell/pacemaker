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


ncfile_path = '/home/nicholat/project/pacemaker/ncfiles/'
temp = iris.load_cube(ncfile_path + 'temp.m48.sfc.4ysl.nc')

# import the ACCESS data using iris
temp_plv = iris.load_cube(ncfile_path + 'temp.m48.plv.4ysl.nc')
temp_plv.coord('p').standard_name = 'air_pressure'

pressure = temp_plv.coord('air_pressure').points

# Define regions
temp_plv.coord('latitude').guess_bounds()
temp_plv.coord('longitude').guess_bounds()

# ---------- define sst/tland----------
lsmask = iris.load_cube(ncfile_path + 'lsmask.nc')[0,0,::]
landmask = ~(ma.make_mask(lsmask.data.copy()) + np.zeros(temp_plv.shape)).astype(bool) # mask sea, show land
seamask = (ma.make_mask(lsmask.data.copy()) + np.zeros(temp_plv.shape)).astype(bool) # mask land, show sea

# sys.exit('exiiiiiiitt')
Tocean = temp_plv.copy()
Tland = temp_plv.copy()
Tocean.data = ma.array(Tocean.data, mask=seamask)
Tland.data = ma.array(Tland.data, mask=landmask)
# --------------
mons = 6
lag = 3
max_i = 35 + lag
max_f = max_i + mons
min_i = 11 + lag
min_f = min_i + mons

plt.clf()
plt.ion()

nlat = 5
latmax = 30
lat_range = np.linspace(0,latmax,nlat).astype('int')
# hold(number of areas, number of lats, pressure levels, max/min)
hold = np.zeros((3,nlat,temp_plv.coord('air_pressure').shape[0],2))
i=0
for lat  in xrange(0,latmax+1,10):
    print(lat)
    NHlati = lat
    NHlatf = lat+10
    SHlati = -lat-10
    SHlatf = -lat

    Tocean_forcNH = Tocean.extract(iris.Constraint(latitude = lambda v: NHlati <= v <= NHlatf, longitude = lambda l: 140 <= l <= 300))
    TO_forc_meanNH = Tocean_forcNH.collapsed(['latitude','longitude'],
                    iris.analysis.MEAN,weights=iris.analysis.cartography.area_weights(Tocean_forcNH))

    Tocean_remNH = Tocean.extract(iris.Constraint(latitude = lambda v: NHlati <= v <= NHlatf,
                    longitude = lambda l: (300 <= l <= 360) or (0 <= l <= 130)))
    TO_rem_meanNH = Tocean_remNH.collapsed(['latitude','longitude'],
                    iris.analysis.MEAN,weights=iris.analysis.cartography.area_weights(Tocean_remNH))

    Tland_latNH = Tland.extract(iris.Constraint(latitude = lambda v: NHlati <= v <= NHlatf))
    TL_meanNH = Tland_latNH.collapsed(['latitude','longitude'],
                    iris.analysis.MEAN,weights=iris.analysis.cartography.area_weights(Tland_latNH))

    Tocean_forcSH = Tocean.extract(iris.Constraint(latitude = lambda v: SHlati <= v <= SHlatf, longitude = lambda l: 140 <= l <= 300))
    TO_forc_meanSH = Tocean_forcSH.collapsed(['latitude','longitude'],
                    iris.analysis.MEAN,weights=iris.analysis.cartography.area_weights(Tocean_forcSH))

    Tocean_remSH = Tocean.extract(iris.Constraint(latitude = lambda v: SHlati <= v <= SHlatf,
                    longitude = lambda l: (300 <= l <= 360) or (0 <= l <= 130)))
    TO_rem_meanSH = Tocean_remSH.collapsed(['latitude','longitude'],
                    iris.analysis.MEAN,weights=iris.analysis.cartography.area_weights(Tocean_remSH))

    Tland_latSH = Tland.extract(iris.Constraint(latitude = lambda v: SHlati <= v <= SHlatf))
    TL_meanSH = Tland_latSH.collapsed(['latitude','longitude'],
                    iris.analysis.MEAN,weights=iris.analysis.cartography.area_weights(Tland_latSH))


    hold[0,i,:,0]   = TO_forc_meanNH[max_i:max_f,::].collapsed('t',iris.analysis.MEAN).data
    hold[0,i,:,1]   = TO_forc_meanSH[max_i:max_f,::].collapsed('t',iris.analysis.MEAN).data
#     hold[0,i,:,1]   = TO_forc_meanNH[min_i:min_f,::].collapsed('t',iris.analysis.MEAN).data
    hold[1,i,:,0]   = TO_rem_meanNH[max_i:max_f,::].collapsed('t',iris.analysis.MEAN).data
    hold[1,i,:,1]   = TO_rem_meanSH[max_i:max_f,::].collapsed('t',iris.analysis.MEAN).data
#     hold[1,i,:,1]   = TO_rem_meanNH[min_i:min_f,::].collapsed('t',iris.analysis.MEAN).data
    hold[2,i,:,0]   = TL_meanNH[max_i:max_f,::].collapsed('t',iris.analysis.MEAN).data
    hold[2,i,:,1]   = TL_meanSH[max_i:max_f,::].collapsed('t',iris.analysis.MEAN).data
#     hold[2,i,:,1]   = TL_meanNH[min_i:min_f,::].collapsed('t',iris.analysis.MEAN).data

    plt.plot(hold[1,i,:,0],pressure,'--',linewidth=2,color=mc.jetloop(i,nlat))
    plt.plot(hold[1,i,:,1]+1,pressure,'--',linewidth=2,color=mc.jetloop(i,nlat))
#     plt.plot(hold[0,i,:,0],pressure,'-.',color=mc.jetloop(i,nlat),label="_nolegend_")
#     plt.plot(hold[0,i,:,1]+1,pressure,'-.',color=mc.jetloop(i,nlat),label="_nolegend_")
    plt.plot(hold[2,i,:,0],pressure,linewidth=2,color=mc.jetloop(i,nlat),label=str(lat))
    plt.plot(hold[2,i,:,1]+1,pressure,linewidth=2,color=mc.jetloop(i,nlat),label=str(lat))

    i=i+1
ext1 = plt.plot(0,0,'-',linewidth=1,color='k',label='Tland')
ext2 = plt.plot(0,0,'--',linewidth=1,color='k',label='Tocean remote')
plt.xlim(-1.2,2.5)
plt.ylabel('z [hPa]')
plt.xlabel('Temp')
plt.gca().invert_yaxis()
plt.legend(loc=3)
plt.title('Temp profile, Tropical Ocean, Max/Min forcing') 
plt.savefig('zonal_profiles.png')

# plt.plot(TOf_max.data,pressure,color='r')
# plt.plot(-TOf_min.data,pressure,'--',color='r',label="_nolegend_")
# plt.plot(TOr_max.data,pressure,color='c')
# plt.plot(-TOr_min.data,pressure,'--',color='c',label="_nolegend_")
# plt.plot(TL_max.data,pressure,color='g')
# plt.plot(-TL_min.data,pressure,'--',color='g',label="_nolegend_")
# plt.legend(('TO','TO_forcing','TO_remote','TL'),loc=0,title='Latitude')

