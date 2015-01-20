import numpy as np
import numpy.ma as ma
import iris as iris
import iris.plot as iplt
import iris.quickplot as qplt
from windspharm.iris import VectorWind
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import mycmaps as mc
import scipy.stats as stats
import sys
import troposave as ta
import mapping as mp
import ntiris as nt

ncfile_path = '/home/nicholat/project/pacemaker/ncfiles/'
u = iris.load_cube(ncfile_path + 'u.plev.4ysl.nc')
v = iris.load_cube(ncfile_path + 'v.plev.4ysl.nc')
tsfc = iris.load_cube(ncfile_path + 'temp.sfc.4ysl.nc')

u.coord('Pressure').standard_name = 'air_pressure'
v.coord('Pressure').standard_name = 'air_pressure'
pressure = u.coord('air_pressure').points

try:
    u.coord('t').standard_name = 'time'
    v.coord('t').standard_name = 'time'
    tsfc.coord('t').standard_name = 'time'
except:
    pass
tsfc = nt.remove_seascyc(tsfc)
v = nt.remove_seascyc(v)
u = nt.remove_seascyc(u)

# Define regions
u.coord('latitude').guess_bounds()
u.coord('longitude').guess_bounds()
v.coord('latitude').guess_bounds()
v.coord('longitude').guess_bounds()
tsfc.coord('latitude').guess_bounds()
tsfc.coord('longitude').guess_bounds()
u.standard_name = 'eastward_wind' 
v.standard_name = 'northward_wind' 

# ---------- define sst/tland----------
lsmask = iris.load_cube(ncfile_path + 'lsmask.nc')[0,0,::]
lsmask.coord('latitude').guess_bounds()
lsmask.coord('longitude').guess_bounds()
landmask = ~(ma.make_mask(lsmask.data.copy()) + np.zeros(tsfc.shape)).astype(bool) # mask sea, show land
seamask = (ma.make_mask(lsmask.data.copy()) + np.zeros(tsfc.shape)).astype(bool) # mask land, show sea

Tocean = tsfc.copy()
Tland = tsfc.copy()
Tocean.data = ma.array(Tocean.data, mask=seamask)
Tland.data = ma.array(Tland.data, mask=landmask)
# landmask = ~(ma.make_mask(lsmask.data.copy()) + np.zeros(temp_plv.shape)).astype(bool) # mask sea, show land
# seamask = (ma.make_mask(lsmask.data.copy()) + np.zeros(temp_plv.shape)).astype(bool) # mask land, show sea

# define upper, middle, lower regions of the troposphere for u, v
u_uptrop = ta.iristropave(u,plev_bottom=500,plev_top=200)
u_midtrop = ta.iristropave(u,plev_bottom=850,plev_top=500)
u_lowtrop = ta.iristropave(u,plev_bottom=1000,plev_top=850)
v_uptrop = ta.iristropave(v,plev_bottom=500,plev_top=200)
v_midtrop = ta.iristropave(v,plev_bottom=850,plev_top=500)
v_lowtrop = ta.iristropave(v,plev_bottom=1000,plev_top=850)
# 
# Define Australian region and generate timeseries of mean Aus sfc temps
temp_sfc_anom = Tland - Tland.collapsed('time',iris.analysis.MEAN)
Aus = iris.Constraint(longitude=lambda l: (109 <= l <= 160), latitude = lambda l: (-42 <= l <= -11))
Tsfc_Aus = temp_sfc_anom.extract(Aus)
grid_areas_Aus = iris.analysis.cartography.area_weights(Tsfc_Aus)

Tsfc_Aus_mean = Tsfc_Aus.collapsed(['latitude', 'longitude'],
                           iris.analysis.MEAN,
                           weights=grid_areas_Aus)


regu_low, coru_low,  pvalu_low = nt.linregts(u_lowtrop,Tsfc_Aus_mean[:,0],'u_lowtrop','tsfc_Aus')
regu_mid, coru_mid,  pvalu_mid = nt.linregts(u_midtrop,Tsfc_Aus_mean[:,0],'u_midtrop','tsfc_Aus')
regu_up,  coru_up,  pvalu_up  = nt.linregts(u_uptrop,Tsfc_Aus_mean[:,0],'u_uptrop','tsfc_Aus')
regv_low, corv_low,  pvalv_low = nt.linregts(v_lowtrop,Tsfc_Aus_mean[:,0],'v_lowtrop','tsfc_Aus')
regv_mid, corv_mid,  pvalv_mid = nt.linregts(v_midtrop,Tsfc_Aus_mean[:,0],'v_midtrop','tsfc_Aus')
regv_up,  corv_up,  pvalv_up  = nt.linregts(v_uptrop,Tsfc_Aus_mean[:,0],'v_uptrop','tsfc_Aus')

qplt.pcmeshclf(regu_low,vmin=-0.2,vmax=0.2,cmap=mc.jetwhite())
plt.title('Regression between Tsfc Aus, u lower tropos')
plt.savefig('./figures/aus_ulow_reg.pdf')
qplt.pcmeshclf(coru_low,vmin=-0.7,vmax=0.7,cmap=mc.jetwhite())
plt.title('Correlation between Tsfc Aus, u lower tropos')
plt.savefig('./figures/aus_ulow_cor.pdf')
qplt.pcmeshclf(regu_mid,vmin=-0.2,vmax=0.2,cmap=mc.jetwhite())
plt.title('Regression between Tsfc Aus, u mider tropos')
plt.savefig('./figures/aus_umid_reg.pdf')
qplt.pcmeshclf(coru_mid,vmin=-0.7,vmax=0.7,cmap=mc.jetwhite())
plt.title('Correlation between Tsfc Aus, u mider tropos')
plt.savefig('./figures/aus_umid_cor.pdf')
qplt.pcmeshclf(regu_up,vmin=-0.2,vmax=0.2,cmap=mc.jetwhite())
plt.title('Regression between Tsfc Aus, u uper tropos')
plt.savefig('./figures/aus_uup_reg.pdf')
qplt.pcmeshclf(coru_up,vmin=-0.7,vmax=0.7,cmap=mc.jetwhite())
plt.title('Correlation between Tsfc Aus, u uper tropos')
plt.savefig('./figures/aus_uup_cor.pdf')

qplt.pcmeshclf(regv_low,vmin=-0.2,vmax=0.2,cmap=mc.jetwhite())
plt.title('Regression between Tsfc Aus, v lower tropos')
plt.savefig('./figures/aus_vlow_reg.pdf')
qplt.pcmeshclf(corv_low,vmin=-0.7,vmax=0.7,cmap=mc.jetwhite())
plt.title('Correlation between Tsfc Aus, v lower tropos')
plt.savefig('./figures/aus_vlow_cor.pdf')
qplt.pcmeshclf(regv_mid,vmin=-0.2,vmax=0.2,cmap=mc.jetwhite())
plt.title('Regression between Tsfc Aus, v mider tropos')
plt.savefig('./figures/aus_vmid_reg.pdf')
qplt.pcmeshclf(corv_mid,vmin=-0.7,vmax=0.7,cmap=mc.jetwhite())
plt.title('Correlation between Tsfc Aus, v mider tropos')
plt.savefig('./figures/aus_vmid_cor.pdf')
qplt.pcmeshclf(regv_up,vmin=-0.2,vmax=0.2,cmap=mc.jetwhite())
plt.title('Regression between Tsfc Aus, v uper tropos')
plt.savefig('./figures/aus_vup_reg.pdf')
qplt.pcmeshclf(corv_up,vmin=-0.7,vmax=0.7,cmap=mc.jetwhite())
plt.title('Correlation between Tsfc Aus, v uper tropos')
plt.savefig('./figures/aus_vup_cor.pdf')

