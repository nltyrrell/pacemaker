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
import pickle

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


def regions(cube,clim=False,maxf=True):
    mons = 4
    lag = 1
    max_i = 35 + lag
    max_f = max_i + mons
    min_i = 11 + lag
    min_f = min_i + mons
    if maxf:
        forc_i = max_i
        forc_f = max_f
    else:
        forc_i = min_i
        forc_f = min_f


    if clim:
        cube_max = cube.collapsed('time',iris.analysis.MEAN)
    else:
        cube_max = cube[forc_i:forc_f,:,::].collapsed('time',iris.analysis.MEAN)
    # cube_anom = cube-cube_mean

    cube_max.coord('latitude').guess_bounds()
    cube_max.coord('longitude').guess_bounds()

    print cube_max.shape
    # ---------- define sst/tland----------
    lsmask = iris.load_cube(ncfile_path + 'lsmask.nc')[0,0,::]
    landmask = ~(ma.make_mask(lsmask.data.copy()) + np.zeros(cube_max.shape)).astype(bool) # mask sea, show land
    # mask out the enso region: masked = 1, not masked = 0
    # mask from x = 43:73, y = 32:40
    enso_mask = np.zeros(cube_max.shape)
    enso_mask[0,30:42,43:75] = 1
    seamask = (ma.make_mask(lsmask.data.copy()) + np.zeros(cube_max.shape) + enso_mask).astype(bool) # mask land, show sea

    Cocean = cube_max.copy()
    Cland = cube_max.copy()
    Cocean.data = ma.array(Cocean.data, mask=seamask)
    Cland.data = ma.array(Cland.data, mask=landmask)
    # --------------
    land_cube = Cland
    ocean_cube = Cocean

    oclat10, oclat10_mean = regmean(ocean_cube,loni=0,lonf=360,lati=-10,latf=10)
    oclatp20, oclatp20_mean = regmean(ocean_cube,loni=0,lonf=360,lati=10,latf=20)
    oclatm20, oclatm20_mean = regmean(ocean_cube,loni=0,lonf=360,lati=-20,latf=-10)
    oclatp30, oclatp30_mean = regmean(ocean_cube,loni=0,lonf=360,lati=20,latf=30)
    oclatm30, oclatm30_mean = regmean(ocean_cube,loni=0,lonf=360,lati=-30,latf=-20)

    lat10, lat10_mean = regmean(land_cube,loni=0,lonf=360,lati=-10,latf=10)
    latp20, latp20_mean = regmean(land_cube,loni=0,lonf=360,lati=10,latf=20)
    latm20, latm20_mean = regmean(land_cube,loni=0,lonf=360,lati=-20,latf=-10)
    latp30, latp30_mean = regmean(land_cube,loni=0,lonf=360,lati=20,latf=30)
    latm30, latm30_mean = regmean(land_cube,loni=0,lonf=360,lati=-30,latf=-20)
#     lat20mean = latp20

    return [oclat10_mean, lat10_mean, oclatp20_mean, latp20_mean, oclatm20_mean, latm20_mean, oclatp30_mean, latp30_mean, oclatm30_mean, latm30_mean]
    
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
lhf = iris.load_cube(ncfile_path + 'lhf.sfc.4ysl.m48.nc')
lhf_reg = regions(lhf)
shf = iris.load_cube(ncfile_path + 'shf.sfc.4ysl.m48.nc')
shf_reg = regions(shf)

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
regarr= np.zeros((12,np.shape(temp_reg)[0]))
for n in xrange(np.shape(temp_reg)[0]):
    print n
    regarr[0,n] = stemp_reg[n].data[0]
    regarr[1,n] = temp_700_reg[n].data
    regarr[2,n] = temp_300_reg[n].data
    regarr[3,n] = rhum_700_reg[n].data
    regarr[4,n] = rhum_300_reg[n].data
    regarr[5,n] = cld_high_reg[n].data
    regarr[6,n] = cld_low_reg[n].data
    regarr[7,n] = precip_reg[n].data[0]
    regarr[8,n] = dlwr_reg[n].data[0]
    regarr[9,n] = dswr_reg[n].data[0]
    regarr[10,n] = lhf_reg[n].data
    regarr[11,n] = shf_reg[n].data
#     regarr[12,n] = u_high_reg[n].data
#     regarr[13,n] = u_low_reg[n].data
#     regarr[14,n] = v_high_reg[n].data
#     regarr[15,n] = v_low_reg[n].data

with open('./pickles/regarr.pickle','wb') as f:
	pickle.dump(regarr,f)

for n in xrange(regarr.shape[0]):
    regarr[n,:] = regarr[n,:]/(regarr[n,:].std())

var = np.array(['Tsfc','T700hPa','T300hpa','RH700hPa','RH300hPa','Cld High','Cld Low','Precip','DLWR','DSWR','lhf','shf']) #,'u_high','u_low','v_high','v_low'])
regs  = np.array(['oclat10', 'lat10', 'oclatN20', 'latN20', 'oclatS20', 'latS20', 'oclatN30', 'latN30', 'oclatS30', 'latS30'])

plt.close('all')
fig, axes = plt.subplots(nrows=1) #,ncols=2)
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
ppl.pcolormesh(fig, axes, regarr, ax_colorbar=cbar_ax, yticklabels=var, xticklabels=regs,vmin=-3.5,vmax=3.5)
axes.set_xlabel('Regions')
axes.set_ylabel('variables')
axes.set_title('Response to Max forcing')
plt.show()
fig.set_size_inches(12,5)
plt.savefig('./figures/lat_response_manyvar_max.pdf')



