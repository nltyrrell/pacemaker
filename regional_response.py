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
import prettyplotlib as ppl

"""
Define many regions, based on surface temperature response (or moisture?)
get the area mean value in these regions for various variables
plot... in some useful mannner tbc
"""

def regmean(cube,loni,lonf,lati,latf):
    """ Define a region and get the area weighted mean
    Input:  cube, lon_i, lon_f, lat_i, lat_f
    Output: cube_reg, cube_regmean
    """
    region = iris.Constraint(longitude=lambda l: (loni <= l <= lonf), latitude = lambda l: (lati <= l <= latf))
    cube_reg = cube.extract(region)
    grid_areas = iris.analysis.cartography.area_weights(cube_reg)
    cube_regmean = cube_reg.collapsed(['latitude', 'longitude'],
                           iris.analysis.MEAN,
                           weights=grid_areas)
#     print 'mean sfc '+str(cube_regmean[0].data)
#     print 'mean p=10 '+str(cube_regmean[10].data)
    return cube_reg, cube_regmean


def regions(cube,clim=False):
    if clim:
        cube_max = cube.collapsed('time',iris.analysis.MEAN)
    else:
        cube_max = cube[35:42,:,::].collapsed('time',iris.analysis.MEAN)
    # cube_anom = cube-cube_mean

    cube_max.coord('latitude').guess_bounds()
    cube_max.coord('longitude').guess_bounds()

    print cube_max.shape
    # ---------- define sst/tland----------
    lsmask = iris.load_cube(ncfile_path + 'lsmask.nc')[0,0,::]
    landmask = ~(ma.make_mask(lsmask.data.copy()) + np.zeros(cube_max.shape)).astype(bool) # mask sea, show land
    seamask = (ma.make_mask(lsmask.data.copy()) + np.zeros(cube_max.shape)).astype(bool) # mask land, show sea

    Cocean = cube_max.copy()
    Cland = cube_max.copy()
    Cocean.data = ma.array(Cocean.data, mask=seamask)
    Cland.data = ma.array(Cland.data, mask=landmask)
    # --------------
    land_cube = Cland

    India, India_mean = regmean(land_cube,loni=60,lonf=75,lati=0,latf=30)
    MC, MC_mean = regmean(land_cube,loni=90,lonf=140,lati=-10,latf=10)
    TropSthAm, TropSthAm_mean = regmean(land_cube,loni=290,lonf=315,lati=-23,latf=0)
    SthSthAm, SthSthAm_mean = regmean(land_cube,loni=270,lonf=315,lati=-60,latf=-24)
    NthWestAfr, NthWestAfr_mean = regmean(land_cube,loni=-15,lonf=15,lati=10,latf=30)
    NthEastAfr, NthEastAfr_mean = regmean(land_cube,loni=15,lonf=50,lati=10,latf=30)
    TropAfr, TropAfr_mean = regmean(land_cube,loni=12,lonf=40,lati=-15,latf=5)
    Aus, Aus_mean = regmean(land_cube,loni=120,lonf=140,lati=-30,latf=-17)

    return [India_mean , MC_mean , TropSthAm_mean , SthSthAm_mean , NthWestAfr_mean , NthEastAfr_mean , TropAfr_mean , Aus_mean]
    
ncfile_path = '/home/nicholat/project/pacemaker/ncfiles/'

temp_plv = iris.load_cube(ncfile_path + 'temp.plv.4ysl.m48.nc')
temp_reg = regions(temp_plv)
temp_sfc = iris.load_cube(ncfile_path + 'temp.sfc.4ysl.m48.nc')
stemp_reg = regions(temp_sfc)
rhum_plv = iris.load_cube(ncfile_path + 'rhum.plv.4ysl.m48.nc')
rhum_reg = regions(rhum_plv)
smc = iris.load_cube(ncfile_path + 'smc.sfc.4ysl.m48.nc')
smc_clim_reg = regions(smc,clim=True)
smc_reg = regions(smc)
dlwr = iris.load_cube(ncfile_path + 'dlwr.sfc.4ysl.m48.nc')
dlwr_reg = regions(dlwr)
dswr = iris.load_cube(ncfile_path + 'dswr.sfc.4ysl.m48.nc')
dswr_reg = regions(dswr)
precip = iris.load_cube(ncfile_path + 'precip.4ysl.m48.nc')
precip_reg = regions(precip)
cld = iris.load_cube(ncfile_path + 'cld.thlev.4ysl.m48.nc')
cld_reg = regions(cld)
u_thlv = iris.load_cube(ncfile_path + 'u.thlev.4ysl.fix.m48.nc')
u_reg = regions(u_thlv)
v_thlv = iris.load_cube(ncfile_path + 'v.thlev.4ysl.m48.nc')
v_thlv = v_thlv.regrid(u_thlv,iris.analysis.Linear())
v_reg = regions(v_thlv)

high = iris.Constraint(atmosphere_hybrid_height_coordinate = lambda h: 5000 <= h <= 15000)
low = iris.Constraint(atmosphere_hybrid_height_coordinate = lambda h: 0 <= h <= 5000)
full = iris.Constraint(atmosphere_hybrid_height_coordinate = lambda h: 0 <= h <= 10000)
cld_high_reg = {} 
cld_low_reg = {} 
for n, i in enumerate(cld_reg):
    i.coord('Hybrid height').standard_name = 'atmosphere_hybrid_height_coordinate'
    cld_high_reg[n] = i.extract(high).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)
    cld_low_reg[n] = i.extract(low).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)

temp_700_reg = {} 
temp_300_reg = {} 
for n, i in enumerate(temp_reg):
    temp_700_reg[n] = i.extract(iris.Constraint(p=700))
    temp_300_reg[n] = i.extract(iris.Constraint(p=300))

u_high_reg = {} 
u_low_reg = {} 
for n, i in enumerate(u_reg):
#     i.coord('Hybrid height').standard_name = 'atmosphere_hybrid_height_coordinate'
    u_high_reg[n] = i.extract(high).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)
    u_low_reg[n] = i.extract(low).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)

v_high_reg = {} 
v_low_reg = {} 
for n, i in enumerate(v_reg):
    i.coord('Hybrid height').standard_name = 'atmosphere_hybrid_height_coordinate'
    v_high_reg[n] = i.extract(high).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)
    v_low_reg[n] = i.extract(low).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)

rhum_700_reg = {} 
rhum_300_reg = {} 
for n, i in enumerate(rhum_reg):
    rhum_700_reg[n] = i.extract(iris.Constraint(p=700))
    rhum_300_reg[n] = i.extract(iris.Constraint(p=300))

# Make array of variables var x reg
regarr= np.zeros((16,np.shape(temp_reg)[0]))
for n in xrange(np.shape(temp_reg)[0]):
    print n
    regarr[0,n] = stemp_reg[n].data[0]
    regarr[1,n] = smc_clim_reg[n].data[0]
    regarr[2,n] = temp_700_reg[n].data
    regarr[3,n] = temp_300_reg[n].data
    regarr[4,n] = rhum_700_reg[n].data
    regarr[5,n] = rhum_300_reg[n].data
    regarr[6,n] = dlwr_reg[n].data[0]
    regarr[7,n] = dswr_reg[n].data[0]
    regarr[8,n] = smc_reg[n].data[0]
    regarr[9,n] = cld_high_reg[n].data
    regarr[10,n] = cld_low_reg[n].data
    regarr[11,n] = precip_reg[n].data[0]
    regarr[12,n] = u_high_reg[n].data
    regarr[13,n] = u_low_reg[n].data
    regarr[14,n] = v_high_reg[n].data
    regarr[15,n] = v_low_reg[n].data

for n in xrange(regarr.shape[0]):
    regarr[n,:] = regarr[n,:]/(regarr[n,:].std())

var = np.array(['Tsfc','smc clim','T700hPa','T300hpa','RH700hPa','RH300hPa','DLWR','DSWR','smc','Cld High','Cld Low','Precip','u_high','u_low','v_high','v_low'])
regs =np.array(['India_mean','MC_mean','TropSthAm_mean','SthSthAm_mean','NthWestAfr_mean','NthEastAfr_mean','TropAfr_mean','Aus_mean'])

plt.close('all')
fig, axes = plt.subplots(nrows=1) #,ncols=2)
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
ppl.pcolormesh(fig, axes, regarr, ax_colorbar=cbar_ax, yticklabels=var, xticklabels=regs,vmin=-2.6,vmax=2.6)
axes.set_xlabel('Regions')
axes.set_ylabel('variables')
axes.set_title('Response to Max forcing')
plt.show()



