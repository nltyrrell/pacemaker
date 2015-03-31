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
var = 'sst'# 'cld'# 'smc'# 'uv'# cloud, smc, or u&v
mm = 'max' # max or min forcing

ncfile_path = '/home/nicholat/project/pacemaker/ncfiles/'
mons = 6
lag = 0
max_i = 35 + lag
max_f = max_i + mons
min_i = 11 + lag
min_f = min_i + mons

if var == 'sst':
    
    sst = iris.load_cube(ncfile_path + 'temp.sfc.4ysl.m48.nc')

    # Define regions
    sst.coord('latitude').guess_bounds()
    sst.coord('longitude').guess_bounds()

    if mm == 'max':
        sst_max            = sst[max_i:max_f,0,::].collapsed('time',iris.analysis.MEAN)
    if mm == 'min':
        sst_max            = sst[min_i:min_f,0,::].collapsed('time',iris.analysis.MEAN)

#     sst_max.units = iris.unit.Unit('K')
    sst_max.long_name = 'sst'
    iris.save(sst_max,'./ncfiles/sst_companom_max.nc')


if var == 'cld':
    # import the ACCESS data using iris
    cld = iris.load_cube(ncfile_path + 'cld.thlev.4ysl.m48.nc')
    cloud_clim = iris.load_cube(ncfile_path + 'cld.thlev.4ysl.nc')
    cld.coord('Hybrid height').standard_name = 'atmosphere_hybrid_height_coordinate'
    cloud_clim.coord('Hybrid height').standard_name = 'atmosphere_hybrid_height_coordinate'

    # Define regions
    cld.coord('latitude').guess_bounds()
    cld.coord('longitude').guess_bounds()
    cloud_clim.coord('latitude').guess_bounds()
    cloud_clim.coord('longitude').guess_bounds()

    high = iris.Constraint(atmosphere_hybrid_height_coordinate = lambda h: 5000 <= h <= 15000)
    low = iris.Constraint(atmosphere_hybrid_height_coordinate = lambda h: 0 <= h <= 5000)
    full = iris.Constraint(atmosphere_hybrid_height_coordinate = lambda h: 0 <= h <= 10000)

    if mm == 'max':
        cld_max = cld[max_i:max_f,::].collapsed('time',iris.analysis.MEAN)
    if mm == 'min':
        cld_max = cld[min_i:min_f,::].collapsed('time',iris.analysis.MEAN)

    cld_high_max = cld_max.extract(high).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.SUM)
    cld_low_max = cld_max.extract(low).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.SUM)
    cld_full_max = cld_max.extract(full).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.SUM)

    cld_high_max.units = iris.unit.Unit(1)
    cld_low_max.units = iris.unit.Unit(1)
    cld_full_max.units = iris.unit.Unit(1)
    cld_high_max.long_name = 'cloud cover'
    cld_low_max.long_name = 'cloud cover'
    cld_full_max.long_name = 'cloud cover'

    iris.save(cld_full_max,'./ncfiles/cld_full_companom_max.nc')
    iris.save(cld_low_max,'./ncfiles/cld_low_companom_max.nc')
    iris.save(cld_high_max,'./ncfiles/cld_high_companom_max.nc')



if var == 'smc':
    
    smc = iris.load_cube(ncfile_path + 'smc.sfc.4ysl.m48.nc')
    smc_clim = iris.load_cube(ncfile_path + 'smc.sfc.4ysl.nc')

    # Define regions
    smc.coord('latitude').guess_bounds()
    smc.coord('longitude').guess_bounds()

    if mm == 'max':
        smc_max            = smc[max_i:max_f,0,::].collapsed('time',iris.analysis.MEAN)
    if mm == 'min':
        smc_max            = smc[min_i:min_f,0,::].collapsed('time',iris.analysis.MEAN)

    smc_max.units = iris.unit.Unit(1)
    smc_max.long_name = 'soil moisture'
    smc_max = smc_max/100
    iris.save(smc_max,'./ncfiles/smc_companom_max.nc')


if var == 'uv':
    u = iris.load_cube(ncfile_path + 'u.plev.4ysl.m48.nc')
    v = iris.load_cube(ncfile_path + 'v.plev.4ysl.m48.nc')

    u.coord('Pressure').standard_name = 'air_pressure'
    v.coord('Pressure').standard_name = 'air_pressure'
    pressure = u.coord('air_pressure').points

    # Define regions
    u.coord('latitude').guess_bounds()
    u.coord('longitude').guess_bounds()
    v.coord('latitude').guess_bounds()
    v.coord('longitude').guess_bounds()
    u.standard_name = 'eastward_wind' 
    v.standard_name = 'northward_wind' 

    plev_850 = iris.Constraint(air_pressure = 850)

    u_850 = u.extract(plev_850)
    if mm == 'max':
        u_850_max = u_850[max_i:max_f,::].collapsed('time',iris.analysis.MEAN)
    if mm == 'min':
        u_850_max = u_850[min_i:min_f,::].collapsed('time',iris.analysis.MEAN)
    iris.save(u_850_max,'./ncfiles/u_850_companom_max.nc')

    v_850 = v.extract(plev_850)
    if mm == 'max':
        v_850_max = v_850[max_i:max_f,::].collapsed('time',iris.analysis.MEAN)
    if mm == 'min':
        v_850_max = v_850[min_i:min_f,::].collapsed('time',iris.analysis.MEAN)
    iris.save(v_850_max,'./ncfiles/v_850_companom_max.nc')



