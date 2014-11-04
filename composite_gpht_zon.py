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
gpht = iris.load_cube(ncfile_path + 'gpht.plv.4ysl.m48.nc')

# import the ACCESS data using iris
# temp_plv = iris.load_cube(ncfile_path + 'theta.plv.4ysl.m48.nc')
# temp_plv.coord('p').standard_name = 'air_pressure'

gpht.coord('latitude').guess_bounds()
gpht.coord('longitude').guess_bounds()
gpht.coord('p').standard_name = 'air_pressure'
gpht.coord('air_pressure').units = 'hPa'
gpht.standard_name = 'geopotential_height_anomaly' 

mons = 6
lag = 0
max_i = 35 + lag
max_f = max_i + mons
min_i = 11 + lag
min_f = min_i + mons

plt.clf()
plt.ion()

lsmask = iris.load_cube(ncfile_path + 'lsmask.nc')[0,0,::]
lsmask.coord('latitude').guess_bounds()
lsmask.coord('longitude').guess_bounds()
# landmask = ~(ma.make_mask(lsmask.data.copy()) + np.zeros(temp_plv.shape)).astype(bool) # mask sea, shogpht land
# seamask = (ma.make_mask(lsmask.data.copy()) + np.zeros(temp_plv.shape)).astype(bool) # mask land, shogpht sea

nlat = 9
latmax = 40
lat_range = np.linspace(0,latmax,nlat).astype('int')
# hold(number of areas, number of lats, pressure levels, max/min)
hold = np.zeros((3,nlat,gpht.coord('air_pressure').shape[0],2))
i=0
for lat  in xrange(0,latmax+1,5):
    print(lat)
    NHlati = lat
    NHlatf = lat+5
    SHlati = -lat-5
    SHlatf = -lat

    gphtNH= gpht.extract(iris.Constraint(latitude = lambda v: NHlati <= v <= NHlatf))
    WNH = gphtNH.collapsed(['latitude'],
                    iris.analysis.MEAN,weights=iris.analysis.cartography.area_weights(gphtNH))
    lsmNH= lsmask.extract(iris.Constraint(latitude = lambda v: NHlati <= v <= NHlatf))
    LSMNH = lsmNH.collapsed(['latitude'],
                    iris.analysis.MEAN,weights=iris.analysis.cartography.area_weights(lsmNH))

    gphtNH_max = WNH[max_i:max_f,::].collapsed('time',iris.analysis.MEAN)
    gphtNH_min = WNH[min_i:min_f,::].collapsed('time',iris.analysis.MEAN)

    gphtSH= gpht.extract(iris.Constraint(latitude = lambda v: SHlati <= v <= SHlatf))
    WSH = gphtSH.collapsed(['latitude'],
                    iris.analysis.MEAN,weights=iris.analysis.cartography.area_weights(gphtSH))

    lsmSH= lsmask.extract(iris.Constraint(latitude = lambda v: SHlati <= v <= SHlatf))
    LSMSH = lsmSH.collapsed(['latitude'],
                    iris.analysis.MEAN,weights=iris.analysis.cartography.area_weights(lsmSH))
    gphtSH_max = WSH[max_i:max_f,::].collapsed('time',iris.analysis.MEAN)
    gphtSH_min = WSH[min_i:min_f,::].collapsed('time',iris.analysis.MEAN)

    gmm = 50
    plt.clf()
    qplt.contourf(gphtNH_max,levels=np.linspace(-gmm,gmm,51),cmap=plt.cm.seismic,extend='both')
    qplt.plot(-20*LSMNH+1000,linewidth=2,color='k')
    latstr = str(gphtNH_max.coord('latitude').points[0])
    plt.title(gphtNH_max.standard_name+' '+latstr+' max forcing')
    plt.xticks(range(0,420,60))
    plt.savefig('./figures/gpht_lat'+latstr+'.png')

    plt.clf()
    qplt.contourf(gphtSH_max,levels=np.linspace(-gmm,gmm,51),cmap=plt.cm.seismic,extend='both')
    latstr = str(gphtSH_max.coord('latitude').points[0])
    qplt.plot(-20*LSMSH+1000,linewidth=2,color='k')
    plt.title(gphtSH_max.standard_name+', lat:'+latstr+', max forcing')
    plt.xticks(range(0,420,60))
    plt.savefig('./figures/gpht_lat'+latstr+'.png')

    gmm = 50
    plt.clf()
    qplt.contourf(gphtNH_min,levels=np.linspace(-gmm,gmm,51),cmap=plt.cm.seismic_r,extend='both')
    qplt.plot(-20*LSMNH+1000,linewidth=2,color='k')
    latstr = str(gphtNH_min.coord('latitude').points[0])
    plt.title(gphtNH_min.standard_name+' '+latstr+' min forcing')
    plt.xticks(range(0,420,60))
    plt.savefig('./figures/gpht_lat'+latstr+'_min.png')

    plt.clf()
    qplt.contourf(gphtSH_min,levels=np.linspace(-gmm,gmm,51),cmap=plt.cm.seismic_r,extend='both')
    latstr = str(gphtSH_min.coord('latitude').points[0])
    qplt.plot(-20*LSMSH+1000,linewidth=2,color='k')
    plt.title(gphtSH_min.standard_name+', lat:'+latstr+', min forcing')
    plt.xticks(range(0,420,60))
    plt.savefig('./figures/gpht_lat'+latstr+'_min.png')

#     hold[0,i,:,0]   = WNH[max_i:max_f,::].collapsed('time',iris.analysis.MEAN).data
#     hold[0,i,:,1]   = TO_forc_meanSH[max_i:max_f,::].collapsed('time',iris.analysis.MEAN).data
#     hold[0,i,:,1]   = TO_forc_meanNH[min_i:min_f,::].collapsed('time',iris.analysis.MEAN).data

#     plt.plot((hold[1,i,:,0]+hold[1,i,:,1])/2,pressure,'--',linewidth=2,color=mc.jetloop(i,nlat))
#     plt.plot((hold[2,i,:,0]+hold[2,i,:,1])/2,pressure,linewidth=2,color=mc.jetloop(i,nlat),label=str(lat))
#     plt.plot(hold[2,i,:,1]+1,pressure,linewidth=2,color=mc.jetloop(i,nlat),label=str(lat))

    i=i+1



# ext1 = plt.plot(0,0,'-',linewidth=1,color='k',label='Tland')
# ext2 = plt.plot(0,0,'--',linewidth=1,color='k',label='Tocean remote')
# plt.ylim(500,1000)
# plt.xlim(-0.0,1.8)
# plt.ylabel('z [hPa]')
# plt.xlabel('Temp')
# plt.gca().invert_yaxis()
# # plt.legend(loc=3)
# plt.title('Zonal Theta anom, Max forcing, '+str(lag)+' to '+str(lag+mons)+' months from peak forcing') 
# plt.savefig('./figures/theta_zonal_prof_zoom.png')

# plt.plot(TOf_max.data,pressure,color='r')
# plt.plot(-TOf_min.data,pressure,'--',color='r',label="_nolegend_")
# plt.plot(TOr_max.data,pressure,color='c')
# plt.plot(-TOr_min.data,pressure,'--',color='c',label="_nolegend_")
# plt.plot(TL_max.data,pressure,color='g')
# plt.plot(-TL_min.data,pressure,'--',color='g',label="_nolegend_")
# plt.legend(('TO','TO_forcing','TO_remote','TL'),loc=0,title='Latitude')

