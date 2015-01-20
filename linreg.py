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


def linregcube(cube1,cube2,name1,name2,ncfile_path='/home/nicholat/project/mit_tcm/access_runs/ncfiles/'):
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

def linregts(cube1,cube2,name1='name1',name2='name2',ncfile_path='/home/nicholat/project/pacemaker/ncfiles/'):
    """
    Calculate the regression between a 2d map and a 1d timeseries
    i.e. regress v onto T_sfc_Aus
    Input: cube_map, cube_ts (time coord only), name1, name2
    Output: cube of regression map
    """

    linreg_map = np.zeros(cube1.shape[1:3])
    cor_map = np.zeros(cube1.shape[1:3])
    #     cor_map = np.zeros(copy_cube.shape)
    #     print('Linreg/Cor map for '+name1+' and '+name2)

    regress_ts = cube2[:,0].data
    for nlat, lat in enumerate(cube1.coord('latitude')):
        for nlon, lon in enumerate(cube1.coord('longitude')):
            # get the sfc temp timeseries at each lat, lon.
            var1 = cube1.extract(iris.Constraint(latitude=lat.points[0],longitude=lon.points[0])).data
    #             var2 = cube2.extract(iris.Constraint(latitude=lat.points[0],longitude=lon.points[0])).data

            linreg = stats.linregress(var1,regress_ts)
            linreg_map[nlat,nlon] = linreg[0]
            cor_map[nlat,nlon] = linreg[2]

        

    reg_cube = cube1.copy()
    reg_cube.data[:] = linreg_map
    reg_cube.long_name = 'Lin Regression '+ name1 +' '+ name2
    reg_cube.units = 'no_unit'
    reg_cube.attributes['title'] = 'Lin Regression '+ name1 +' '+ name2
    reg_cube.attributes['name'] = 'reg'
    reg_cube.remove_coord('surface')
    reg_cube.remove_coord('time')
    iris.save(reg_cube,ncfile_path+'lreg.4ysl.'+name1+'.'+name2+'.nc')
    cor_cube = cube1.copy()
    cor_cube.data[:] = cor_map
    cor_cube.long_name = 'Correlation '+ name1 +' '+ name2
    cor_cube.units = 'no_unit'
    cor_cube.attributes['title'] = 'Correlation '+ name1 +' '+ name2
    cor_cube.attributes['name'] = 'r_val'
    cor_cube.remove_coord('surface')
    cor_cube.remove_coord('time')
    iris.save(cor_cube,ncfile_path+'cor.4ysl.'+name1+'.'+name2+'.nc')

    return reg_cube, cor_cube

# ncfile_path = '/home/nicholat/project/pacemaker/ncfiles/'
# temp = iris.load_cube(ncfile_path + 'temp.sfc.4ysl.ym.nc')
# temp_mean = temp[:,0,::].collapsed('time',iris.analysis.MEAN)
# 
# linregcube(T_sfc,T_tropo,'T_sfc','T_tropo')






