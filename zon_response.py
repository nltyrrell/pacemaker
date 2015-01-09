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

def zonmean(cube):
    """ Define a region and get the area weighted mean
    Input:  cube
    Output: cube_reg, cube_regmean
    """
#     region = iris.Constraint(longitude=lambda l: (loni <= l <= lonf), latitude = lambda l: (lati <= l <= latf))
#     cube_reg = cube.extract(region)
#     grid_areas = iris.analysis.cartography.area_weights(cube_reg)
    cube_zonmean = cube.collapsed(['longitude'],
                           iris.analysis.MEAN)
#                            weights=grid_areas)
#     print 'mean sfc '+str(cube_regmean[0].data)
#     print 'mean p=10 '+str(cube_regmean[10].data)
    return cube_zonmean


def regions(cube,clim=False,maxf=True, maskprint=False):
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
    enso_mask[:,20:56,43:75] = 1
    seamask = (ma.make_mask(lsmask.data.copy()) + np.zeros(cube_max.shape) + enso_mask).astype(bool) # mask land, show sea
#     if maskprint:


    ocean_cube = cube_max.copy()
    land_cube = cube_max.copy()
    ocean_cube.data = ma.array(ocean_cube.data, mask=seamask)
    land_cube.data = ma.array(land_cube.data, mask=landmask)
    # --------------
    oc_zonmean = ocean_cube.collapsed(['longitude'],
                           iris.analysis.MEAN)
    ln_zonmean = land_cube.collapsed(['longitude'],
                           iris.analysis.MEAN)

    return [oc_zonmean, ln_zonmean]
    
ncfile_path = '/home/nicholat/project/pacemaker/ncfiles/'

temp_plv = iris.load_cube(ncfile_path + 'temp.plv.4ysl.m48.nc')
temp_oc, temp_ln = regions(temp_plv)
temp_sfc = iris.load_cube(ncfile_path + 'temp.sfc.4ysl.m48.nc')
stemp_oc, stemp_ln = regions(temp_sfc)
rhum_plv = iris.load_cube(ncfile_path + 'rhum.plv.4ysl.m48.nc')
rhum_oc, rhum_ln = regions(rhum_plv)
smc = iris.load_cube(ncfile_path + 'smc.sfc.4ysl.m48.nc')
smc_clim_oc, smc_clim_ln = regions(smc,clim=True)
smc_oc, smc_ln = regions(smc)
dlwr = iris.load_cube(ncfile_path + 'dlwr.sfc.4ysl.m48.nc')
dlwr_oc, dlwr_ln = regions(dlwr)
dswr = iris.load_cube(ncfile_path + 'dswr.sfc.4ysl.m48.nc')
dswr_oc, dswr_ln = regions(dswr)
precip = iris.load_cube(ncfile_path + 'precip.4ysl.m48.nc')
precip_oc, precip_ln = regions(precip)
cld = iris.load_cube(ncfile_path + 'cld.thlev.4ysl.m48.nc')
cld_oc, cld_ln = regions(cld)
u_thlv = iris.load_cube(ncfile_path + 'u.thlev.4ysl.fix.m48.nc')
u_oc, u_ln = regions(u_thlv)
v_thlv = iris.load_cube(ncfile_path + 'v.thlev.4ysl.m48.nc')
v_thlv = v_thlv.regrid(u_thlv,iris.analysis.Linear())
v_oc, v_ln = regions(v_thlv)
lhf = iris.load_cube(ncfile_path + 'lhf.sfc.4ysl.m48.nc')
lhf_oc, lhf_ln = regions(lhf)
shf = iris.load_cube(ncfile_path + 'shf.sfc.4ysl.m48.nc')
shf_oc, shf_ln = regions(shf)


high = iris.Constraint(atmosphere_hybrid_height_coordinate = lambda h: 5000 <= h <= 15000)
low = iris.Constraint(atmosphere_hybrid_height_coordinate = lambda h: 0 <= h <= 5000)
full = iris.Constraint(atmosphere_hybrid_height_coordinate = lambda h: 0 <= h <= 10000)

cld_oc.coord('Hybrid height').standard_name = 'atmosphere_hybrid_height_coordinate'
cld_high_oc = cld_oc.extract(high).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)
cld_low_oc = cld_oc.extract(low).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)
cld_ln.coord('Hybrid height').standard_name = 'atmosphere_hybrid_height_coordinate'
cld_high_ln = cld_ln.extract(high).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)
cld_low_ln = cld_ln.extract(low).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)

temp_700_oc = temp_oc.extract(iris.Constraint(p=700))
temp_300_oc = temp_oc.extract(iris.Constraint(p=300))
temp_700_ln = temp_ln.extract(iris.Constraint(p=700))
temp_300_ln = temp_ln.extract(iris.Constraint(p=300))

#     i.coord('Hybrid height').standard_name = 'atmosphere_hybrid_height_coordinate'
u_high_oc = u_oc.extract(high).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)
u_low_oc = u_oc.extract(low).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)
u_high_ln = u_ln.extract(high).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)
u_low_ln = u_ln.extract(low).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)

v_oc.coord('Hybrid height').standard_name = 'atmosphere_hybrid_height_coordinate'
v_high_oc = v_oc.extract(high).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)
v_low_oc = v_oc.extract(low).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)
v_ln.coord('Hybrid height').standard_name = 'atmosphere_hybrid_height_coordinate'
v_high_ln = v_ln.extract(high).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)
v_low_ln = v_ln.extract(low).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)

rhum_700_oc = rhum_oc.extract(iris.Constraint(p=700))
rhum_300_oc = rhum_oc.extract(iris.Constraint(p=300))
rhum_700_ln = rhum_ln.extract(iris.Constraint(p=700))
rhum_300_ln = rhum_ln.extract(iris.Constraint(p=300))

# Make array of variables var x reg
regarr= np.zeros((12,2,np.shape(temp_oc)[-1]))
regarr[0,1,:] = stemp_oc.data[0]
regarr[1,1,:] = temp_700_oc.data
regarr[2,1,:] = temp_300_oc.data
regarr[3,1,:] = rhum_700_oc.data
regarr[4,1,:] = rhum_300_oc.data
regarr[5,1,:] = cld_high_oc.data
regarr[6,1,:] = cld_low_oc.data
regarr[7,1,:] = precip_oc.data[0]*1000
regarr[8,1,:] = dlwr_oc.data[0]
regarr[9,1,:] = dswr_oc.data[0]
regarr[10,1,:] = lhf_oc.data
regarr[11,1,:] = shf_oc.data
regarr[0,0,:] = stemp_ln.data[0]
regarr[1,0,:] = temp_700_ln.data
regarr[2,0,:] = temp_300_ln.data
regarr[3,0,:] = rhum_700_ln.data
regarr[4,0,:] = rhum_300_ln.data
regarr[5,0,:] = cld_high_ln.data
regarr[6,0,:] = cld_low_ln.data
regarr[7,0,:] = precip_ln.data[0]*1000
regarr[8,0,:] = dlwr_ln.data[0]
regarr[9,0,:] = dswr_ln.data[0]
regarr[10,0,:] = lhf_ln.data
regarr[11,0,:] = shf_ln.data
#     regarr[12,n] = u_high_reg[n].data
#     regarr[13,n] = u_low_reg[n].data
#     regarr[14,n] = v_high_reg[n].data
#     regarr[15,n] = v_low_reg[n].data

# with open('./pickles/regarr.pickle','wb') as f:
# 	pickle.dump(regarr,f)

# for n in xrange(regarr.shape[0]):
#     regarr[n,0,:] = regarr[n,0,:]/(regarr[n,0,:].std())
#     regarr[n,1,:] = regarr[n,1,:]/(regarr[n,1,:].std())

var = np.array(['Tsfc','T700hPa','T300hpa','RH700hPa','RH300hPa','Cld High','Cld Low','Precip','DLWR','DSWR','lhf','shf']) #,'u_high','u_low','v_high','v_low'])
regs  = np.array(['oclat10', 'lat10', 'oclatN20', 'latN20', 'oclatS20', 'latS20', 'oclatN30', 'latN30', 'oclatS30', 'latS30'])
lats=np.linspace(-90,90,73)

plt.clf()
# fig, axes = plt.subplots(nrows=1) #,ncols=2)
# fig.subplots_adjust(right=0.8)
# cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])

def jet(i,nloops):
    N_jet = plt.cm.jet.N
    colornum = plt.cm.jet((N_jet*i/(nloops-1)))
    return colornum

for n, i in enumerate(regarr[0:3]):
    plt.plot(lats,i[0,:], label=var[n]+', land',color=jet(n,3),linewidth=2)
    plt.plot(lats,i[1,:],'--', label=var[n]+', ocean',color=jet(n,3),linewidth=2)
plt.title('Response to Max forcing')
plt.xlabel('Latitude')
plt.ylabel('Temp anomaly [K]')
plt.legend()
plt.show()
plt.savefig('./figures/zon_temp.pdf')

plt.clf()

for n, i in enumerate(regarr[3:5]):
    plt.plot(lats,i[0,:], label=var[n+3]+', land',color=jet(n,3),linewidth=2)
    plt.plot(lats,i[1,:],'--', label=var[n+3]+', ocean',color=jet(n,3),linewidth=2)
plt.title('Response to Max forcing')
plt.xlabel('Latitude')
plt.ylabel('RH anomaly')
plt.legend()
plt.show()
plt.savefig('./figures/zon_rh.pdf')

plt.clf()
for n, i in enumerate(regarr[5:8]):
    plt.plot(lats,i[0,:], label=var[n+5]+', land',color=jet(n,3),linewidth=2)
    plt.plot(lats,i[1,:],'--', label=var[n+5]+', ocean',color=jet(n,3),linewidth=2)
plt.title('Response to Max forcing')
plt.xlabel('Latitude')
plt.ylabel('Cloud anomaly')
plt.legend()
plt.show()
plt.savefig('./figures/zon_cldprc.pdf')

plt.clf()
for n, i in enumerate(regarr[8:10]):
    plt.plot(lats,i[0,:], label=var[n+8]+', land',color=jet(n,3),linewidth=2)
    plt.plot(lats,i[1,:],'--', label=var[n+8]+', ocean',color=jet(n,3),linewidth=2)
plt.title('Response to Max forcing')
plt.xlabel('Latitude')
plt.ylabel('Rad anomaly')
plt.legend()
plt.show()
plt.savefig('./figures/zon_rad.pdf')

plt.clf()
for n, i in enumerate(regarr[10:12]):
    plt.plot(lats,i[0,:], label=var[n+10]+', land',color=jet(n,3),linewidth=2)
    plt.plot(lats,i[1,:],'--', label=var[n+10]+', ocean',color=jet(n,3),linewidth=2)
plt.title('Response to Max forcing')
plt.xlabel('Latitude')
plt.ylabel('anomaly')
plt.legend()
plt.show()
plt.savefig('./figures/zon_flux.pdf')




