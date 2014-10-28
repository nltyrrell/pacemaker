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

def linregcube(cube1,cube2,name1,name2,ncfile_path='/home/nicholat/project/pacemaker/ncfiles/',copy_cube=temp_mean):
    linreg_map = np.zeros(copy_cube.shape)
    cor_map = np.zeros(copy_cube.shape)
    print('Linreg/Cor map for '+name1+' and '+name2)
    
    for nlat, lat in enumerate(copy_cube.coord('latitude')):
        for nlon, lon in enumerate(copy_cube.coord('longitude')):
            # get the sfc temp timeseries at each lat, lon.
            var1 = cube1.extract(iris.Constraint(latitude=lat.points[0],longitude=lon.points[0])).data
            var2 = cube2.extract(iris.Constraint(latitude=lat.points[0],longitude=lon.points[0])).data

            linreg = stats.linregress(var1,var2)
            linreg_map[nlat,nlon] = linreg[0]
            cor_map[nlat,nlon] = linreg[2]

    #Copy the soil moisture cube then change everything, clean it up, save phi as an netcdf
    reg_cube = copy_cube.copy()
    reg_cube.data[:] = linreg_map
    reg_cube.long_name = 'Lin Regression '+ name1 +' '+ name2
    reg_cube.units = 'no_unit'
    reg_cube.attributes['title'] = 'Lin Regression '+ name1 +' '+ name2
    reg_cube.attributes['name'] = 'reg'
    reg_cube.remove_coord('surface')
    reg_cube.remove_coord('time')
    iris.save(reg_cube,ncfile_path+'lreg.4ysl.'+name1+'.'+name2+'.nc')
    #Copy the soil moisture cube then change everything, clean it up, save phi as an netcdf
    cor_cube = copy_cube.copy()
    cor_cube.data[:] = cor_map
    cor_cube.long_name = 'Correlation '+ name1 +' '+ name2
    cor_cube.units = 'no_unit'
    cor_cube.attributes['title'] = 'Correlation '+ name1 +' '+ name2
    cor_cube.attributes['name'] = 'r_val'
    cor_cube.remove_coord('surface')
    cor_cube.remove_coord('time')
    iris.save(cor_cube,ncfile_path+'cor.4ysl.'+name1+'.'+name2+'.nc')

    return reg_cube, cor_cube

# import the ACCESS data using iris
temp_plv = iris.load_cube(ncfile_path + 'temp.plv.4ysl.ym.nc')

temp_plv_mean = temp_plv[:,:,::].collapsed('time',iris.analysis.MEAN)
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

for np, p in enumerate(temp_plv.coord('air_pressure')):

        Tplv = Tocean_trop_mean.extract(iris.Constraint(air_pressure=p.points[0])).data
        TOsfc = Tocean_trop_mean.extract(iris.Constraint(air_pressure=1000)).data
        linreg_ln = stats.linregress(TOsfc,Tplv)
        TrTochold[np] = linreg_ln[0]

        Tplv = Tplv_SthAm_armean.extract(iris.Constraint(air_pressure=p.points[0])).data
        Tsfc = Tplv_SthAm_armean.extract(iris.Constraint(air_pressure=1000)).data
        TOsfc = Tocean_trop_mean.extract(iris.Constraint(air_pressure=1000)).data
        linreg_ln = stats.linregress(TOsfc,Tplv)
        SthAmhold[0,np] = linreg_ln[0]
        linreg_ln = stats.linregress(Tsfc,Tplv)
        SthAmhold[1,np] = linreg_ln[0]

        Tplv = Tplv_NthAfr_armean.extract(iris.Constraint(air_pressure=p.points[0])).data
        Tsfc = Tplv_NthAfr_armean.extract(iris.Constraint(air_pressure=1000)).data
        linreg_ln = stats.linregress(TOsfc,Tplv)
        NthAfrhold[0,np] = linreg_ln[0]
        linreg_ln = stats.linregress(Tsfc,Tplv)
        NthAfrhold[1,np] = linreg_ln[0]

        Tplv = Tplv_CntAfr_armean.extract(iris.Constraint(air_pressure=p.points[0])).data
        Tsfc = Tplv_CntAfr_armean.extract(iris.Constraint(air_pressure=1000)).data
        linreg_ln = stats.linregress(TOsfc,Tplv)
        CntAfrhold[0,np] = linreg_ln[0]
        linreg_ln = stats.linregress(Tsfc,Tplv)
        CntAfrhold[1,np] = linreg_ln[0]

        Tplv = Tplv_Aus_armean.extract(iris.Constraint(air_pressure=p.points[0])).data
        Tsfc = Tplv_Aus_armean.extract(iris.Constraint(air_pressure=1000)).data
        linreg_ln = stats.linregress(TOsfc,Tplv)
        Aushold[0,np] = linreg_ln[0]
        linreg_ln = stats.linregress(Tsfc,Tplv)
        Aushold[1,np] = linreg_ln[0]

pressure = temp_plv.coord('air_pressure').points
plt.clf()
plt.plot(NthAfrhold[0],pressure,'g')#linestyle='--',color=wetdry(k,nloops))
plt.plot(CntAfrhold[0],pressure,'b')#linestyle='--',color=wetdry(k,nloops))
plt.plot(SthAmhold[0],pressure,'r')#linestyle='--',color=wetdry(k,nloops))
plt.plot(Aushold[0],pressure,'k')#linestyle='--',color=wetdry(k,nloops))
plt.plot(TrTochold,pressure,'c')#linestyle='--',color=wetdry(k,nloops))
plt.plot(NthAfrhold[1],pressure,'--g')#linestyle='--',color=wetdry(k,nloops))
plt.plot(CntAfrhold[1],pressure,'--b')#linestyle='--',color=wetdry(k,nloops))
plt.plot(SthAmhold[1],pressure,'--r')#linestyle='--',color=wetdry(k,nloops))
plt.plot(Aushold[1],pressure,'--k')#linestyle='--',color=wetdry(k,nloops))
#plt.plot(lnochold,pressure,'b')#color=wetdry(k,nloops))
plt.legend(('NthAfr','Central Africa','Sth America','Australia','Tropical Ocean'),loc=3, )
plt.title('Regression with surface temp and temp profile.') 
plt.xlim(-0.5,3.0)
plt.ylabel('z [hPa]')
plt.xlabel('Regresion coeff.')
plt.gca().invert_yaxis()
# plt.savefig('./figures/regprof_cont.png')




