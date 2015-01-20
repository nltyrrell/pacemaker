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
import ntiris as nt
"""
Define many regions, based on surface temperature response (or moisture?)
get the area mean value in these regions for various variables
plot... in some useful mannner tbc
"""

def regmean(cube,loni,lonf,lati,latf,wrap=False):
    """ Define a region and get the area weighted mean
    Input:  cube, lon_i, lon_f, lat_i, lat_f
    Output: cube_reg, cube_regmean
    kwargs: wrap=True if longitude wraps around 0, then loni should be negative, lonf positive
    """
    if not wrap:
        region = iris.Constraint(longitude=lambda l: (loni <= l <= lonf), latitude = lambda l: (lati <= l <= latf))
    if wrap:
        region = iris.Constraint(longitude=lambda l: (0 <= l <= lonf or (360+loni) <= l <= 360), latitude = lambda l: (lati <= l <= latf))
    cube_region = cube.extract(region)
    grid_areas = iris.analysis.cartography.area_weights(cube_region)
    cube_mean = cube_region.collapsed(['latitude', 'longitude'],
                           iris.analysis.MEAN,
                           weights=grid_areas)
#     print 'mean sfc '+str(cube_regmean[0].data)
#     print 'mean p=10 '+str(cube_regmean[10].data)
    return cube_region, cube_mean


def regions(cube,tsfc_cube,name2='name2'):
    """
    Calculates a regression between variable and tsfc for different regions
    Input: cube, tsfc_cube
    Output: array with all regions
    """
    cube = nt.remove_seascyc(cube)
    tsfc_cube = nt.remove_seascyc(tsfc_cube)

    try:
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()
        tsfc_cube.coord('latitude').guess_bounds()
        tsfc_cube.coord('longitude').guess_bounds()
    except:
        pass

    # ---------- define sst/tland----------
    lsmask = iris.load_cube(ncfile_path + 'lsmask.nc')[0,0,::]
    landmask = ~(ma.make_mask(lsmask.data.copy()) + np.zeros(cube.shape)).astype(bool) # mask sea, show land
    seamask = (ma.make_mask(lsmask.data.copy()) + np.zeros(cube.shape)).astype(bool) # mask land, show sea

    ocean_cube = cube.copy()
    land_cube = cube.copy()
    ocean_cube.data = ma.array(ocean_cube.data, mask=seamask)
    land_cube.data = ma.array(land_cube.data, mask=landmask)
    # --------------
    rorc = 1 # regression or correlation 0 - reg, 1, cor

    print "Calc reg/cor for India"
    India, India_mean = regmean(land_cube,loni=60,lonf=75,lati=0,latf=30)
    name1 = 'India tsfc'
    India_reg = nt.linregts(cube[:,0,::],India_mean[:,0],name1,name2)[rorc]
    print "Calc reg/cor for MC"
    MC, MC_mean = regmean(land_cube,loni=90,lonf=140,lati=-10,latf=10)
    name1 = 'MC tsfc'
    MC_reg = nt.linregts(cube[:,0,::],MC_mean[:,0],name1,name2)[rorc]
    print "Calc reg/cor for TropSthAm"
    TropSthAm, TropSthAm_mean = regmean(land_cube,loni=290,lonf=315,lati=-23,latf=0)
    name1 = 'TropSthAm tsfc'
    TropSthAm_reg = nt.linregts(cube[:,0,::],TropSthAm_mean[:,0],name1,name2)[rorc]
    print "Calc reg/cor for SthSthAm"
    SthSthAm, SthSthAm_mean = regmean(land_cube,loni=270,lonf=315,lati=-60,latf=-24)
    name1 = 'SthSthAm tsfc'
    SthSthAm_reg = nt.linregts(cube[:,0,::],SthSthAm_mean[:,0],name1,name2)[rorc]
    print "Calc reg/cor for NthWestAfr"
    NthWestAfr, NthWestAfr_mean = regmean(land_cube,loni=-15,lonf=15,lati=10,latf=30,wrap=True)
    name1 = 'NthWestAfr tsfc'
    NthWestAfr_reg = nt.linregts(cube[:,0,::],NthWestAfr_mean[:,0],name1,name2)[rorc]
    print "Calc reg/cor for NthEastAfr"
    NthEastAfr, NthEastAfr_mean = regmean(land_cube,loni=15,lonf=50,lati=10,latf=30)
    name1 = 'NthEastAfr tsfc'
    NthEastAfr_reg = nt.linregts(cube[:,0,::],NthEastAfr_mean[:,0],name1,name2)[rorc]
    print "Calc reg/cor for TropAfr"
    TropAfr, TropAfr_mean = regmean(land_cube,loni=12,lonf=40,lati=-15,latf=5)
    name1 = 'TropAfr tsfc'
    TropAfr_reg = nt.linregts(cube[:,0,::],TropAfr_mean[:,0],name1,name2)[rorc]
    print "Calc reg/cor for SthAfr"
    SthAfr, SthAfr_mean = regmean(land_cube,loni=12,lonf=40,lati=-35,latf=-15)
    name1 = 'SthAfr tsfc'
    SthAfr_reg = nt.linregts(cube[:,0,::],SthAfr_mean[:,0],name1,name2)[rorc]
    print "Calc reg/cor for Aus"
    Aus, Aus_mean = regmean(land_cube,loni=120,lonf=140,lati=-30,latf=-17)
    name1 = 'Aus tsfc'
    Aus_reg = nt.linregts(cube[:,0,::],Aus_mean[:,0],name1,name2)[rorc]


    return [India_reg] #, MC_reg, TropSthAm_reg , SthSthAm_reg , NthWestAfr_reg , NthEastAfr_reg, TropAfr_reg, SthAfr_reg, Aus_reg]
    
ncfile_path = '/home/nicholat/project/pacemaker/ncfiles/'

# temp_plv = iris.load_cube(ncfile_path + 'temp.plv.4ysl.m48.nc')
# temp_reg = regions(temp_plv)
temp_sfc = iris.load_cube(ncfile_path + 'temp.sfc.4ysl.nc')
# stemp_reg = regions(temp_sfc,temp_sfc,name2='sfc temp')

rhum_plv = iris.load_cube(ncfile_path + 'rhum.plv.4ysl.nc')
rhum_700 = rhum_plv.extract(iris.Constraint(p=700))
rhum_300 = rhum_plv.extract(iris.Constraint(p=300))
rhum_700_reg = regions(rhum_700,temp_sfc,name2='rhum')
rhum_300_reg = regions(rhum_300,temp_sfc,name2='rhum')

smc = iris.load_cube(ncfile_path + 'smc.sfc.4ysl.nc')
smc_reg = regions(smc,temp_sfc,name2='smc')
dlwr = iris.load_cube(ncfile_path + 'dlwr.sfc.4ysl.nc')
dlwr_reg = regions(dlwr,temp_sfc,name2='dlwr')
dswr = iris.load_cube(ncfile_path + 'dswr.sfc.4ysl.nc')
dswr_reg = regions(dswr,temp_sfc)
precip = iris.load_cube(ncfile_path + 'precip.4ysl.nc')
precip_reg = regions(precip,temp_sfc)
cld = iris.load_cube(ncfile_path + 'cld.thlev.4ysl.nc')
cld_reg = regions(cld,temp_sfc)
u_thlv = iris.load_cube(ncfile_path + 'u.thlev.4ysl.fix.nc')
u_reg = regions(u_thlv,temp_sfc)
v_thlv = iris.load_cube(ncfile_path + 'v.thlev.4ysl.nc')
v_thlv = v_thlv.regrid(u_thlv,iris.analysis.Linear())
v_reg = regions(v_thlv,temp_sfc)
lhf = iris.load_cube(ncfile_path + 'lhf.sfc.4ysl.nc')
lhf_reg = regions(lhf,temp_sfc)
shf = iris.load_cube(ncfile_path + 'shf.sfc.4ysl.nc')
shf_reg = regions(shf,temp_sfc)

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



for n in xrange(regarr.shape[0]):
    regarr[n,:] = regarr[n,:]/(regarr[n,:].std())

var = np.array(['Tsfc','smc clim','T700hPa','T300hpa','RH700hPa','RH300hPa','DLWR','DSWR','smc','Cld High','Cld Low','Precip','lhf','shf']) #,'u_high','u_low','v_high','v_low'])
regs =np.array(['India','MC','TropSthAm','SthSthAm','NthWestAfr','NthEastAfr','TropAfr','SthAfr','Aus','lat10', 'latN20', 'latS20', 'latN30', 'latS30'])

plt.close('all')
fig, axes = plt.subplots(nrows=1) #,ncols=2)
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
ppl.pcolormesh(fig, axes, regarr, ax_colorbar=cbar_ax, yticklabels=var, xticklabels=regs,vmin=-2.6,vmax=2.6)
axes.set_xlabel('Regions')
axes.set_ylabel('variables')
axes.set_title('Response to Max forcing')
plt.show()
fig.set_size_inches(18,5)
plt.savefig('./figures/regional_response_manyvar_max.eps')



