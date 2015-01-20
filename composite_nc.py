import iris
import numpy as np
import iris.coord_categorisation
import iris.plot as iplt
import iris.quickplot as qplt
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import sys



def bc(aux_cube, dim_cube, comparable_coord):
    if type(aux_cube.coord(comparable_coord)) is not iris.coords.AuxCoord:
        raise TypeError
    if type(dim_cube.coord(comparable_coord)) is not iris.coords.AuxCoord:
        raise TypeError
    
    aux_cube_coord = aux_cube.coord(comparable_coord)
    dim_cube_coord = dim_cube.coord(comparable_coord)
    aux_cube_dim, = aux_cube.coord_dims(aux_cube_coord)
    dim_cube_dim, = aux_cube.coord_dims(dim_cube_coord)
    
    s_aux_cube = [slice(None)]*len(aux_cube.shape)
    s_aux_cube[aux_cube_dim] = 0
    s_dim_cube = [slice(None)]*len(dim_cube.shape)
    s_dim_cube[dim_cube_dim] = 0
    a = aux_cube[tuple(s_aux_cube)]
    a.attributes = None
    a.cell_methods = None
    b = dim_cube[tuple(s_dim_cube)]
    b.attributes = None
    b.cell_methods = None
    if not a.is_compatible(b):
        iris.util.describe_diff(a, b)
        raise RuntimeError("Cubes are not compatible")
    
    ind = []
    for p in aux_cube.coord(comparable_coord).points:
        i = np.where(dim_cube.coord(comparable_coord).points == p)
        ind.append(i[0][0])

    s = [slice(None)]*len(dim_cube.shape)
    s[dim_cube_dim] = ind
    new_data = dim_cube.data[tuple(s)]
    new_cube = aux_cube.copy()
    new_cube.data = new_data
    new_cube.history = "%s comparable to %s in terms of %s" % (dim_cube.name(),
                                                               aux_cube.name(),
                                                               comparable_coord)
    
    return new_cube

def remove_seascyc(cube, time_name='time'):
    iris.coord_categorisation.add_month_number(cube, 'time', 'month_number')
    cube_mean = cube[:,:,::].collapsed('time',iris.analysis.MEAN)
    print "cube"
    print cube
    print "cube_mean"
    print cube_mean
    cube_anom = cube-cube_mean
    print "cube_anom"
    print cube_anom
    cube_mon_mean = cube_anom.aggregated_by('month_number', iris.analysis.MEAN)
    print "cube_mon_mean"
    print cube_mon_mean
    seasonal_cycle = bc(cube_anom, cube_mon_mean, 'month_number')
    cube_rsc = cube_anom - seasonal_cycle
    return cube_rsc

def enscyc_ag(cube):
    ens = np.tile(np.linspace(1,48,48),24)
    #trim cube
    print cube.shape
    cube = cube[0:ens.shape[0]]
    print cube.shape
    print cube
    cube.coord('month_number').long_name = '48_months'
    cube.coord('48_months').points = ens

    m48 = cube.aggregated_by('48_months',iris.analysis.MEAN)
    return m48
    

def composite_m48(cube_name, ncfile_path_in='/home/nicholat/project/pacemaker/ncfiles/', ncfile_path_out='/home/nicholat/project/pacemaker/ncfiles/',notanom=False):
    cube = iris.load_cube(ncfile_path_in+cube_name)
    print('Loaded cube: '+cube_name)
    try:
        cube.coord('t').standard_name='time'
    except:
        pass
    else:
        print "t coord changed to time"

    if notanom:
        cube_rsc = cube
        iris.coord_categorisation.add_month_number(cube_rsc, 'time', 'month_number')
    else:
        cube_rsc = remove_seascyc(cube) 
    cube_m48 = enscyc_ag(cube_rsc)
    cube_m48.long_name = cube.long_name
    if notanom:
        new_name = cube_name[:-2]+'m48.abs.nc'
    else:
        new_name = cube_name[:-2]+'m48.nc'
    print new_name
    iris.save(cube_m48,ncfile_path_out+new_name)
    return cube_m48, cube_rsc, cube

# u_cube = composite_m48('u.thlev.4ysl.fix.nc')
lhf = composite_m48('lhf.sfc.4ysl.nc')
shf = composite_m48('shf.sfc.4ysl.nc')

# sys.exit('all done')
temp_cube, temp_rsc, cube = composite_m48('temp.sfc.4ysl.nc')
smc = composite_m48('smc.sfc.4ysl.nc')
dswr = composite_m48('dswr.sfc.4ysl.nc', ncfile_path_in='/home/nicholat/project/access/slabosc/ncfiles/')
dlwr = composite_m48('dlwr.sfc.4ysl.nc', ncfile_path_in='/home/nicholat/project/access/slabosc/ncfiles/')
rain_conv = composite_m48('rain.conv.4ysl.nc')
rain_lgscl = composite_m48('rain.lgscale.4ysl.nc')
precip = composite_m48('precip.4ysl.nc')
smc_cube, smc_rsc, cube = composite_m48('smc.plv.4ysl.nc',notanom=True)
rhum_cube, rhum_rsc, cube = composite_m48('rhum.plv.4ysl.nc',notanom=True)
temp_cube, temp_rsc, cube = composite_m48('temp.plv.4ysl.nc',notanom=True)
gpht_cube = composite_m48('gpht.plv.4ysl.nc')
rh_cube = composite_m48('rhum.plv.4ysl.nc')
p_cube = composite_m48('pres.sfc.4ysl.nc')
v_cube = composite_m48('v.plev.4ysl.nc')
u_cube = composite_m48('u.plev.4ysl.nc')
w_cube = composite_m48('w.thlev.4ysl.fix.nc')
composite_m48('lwflux.clsky.sfc.4ysl.nc')
composite_m48('lhf.sfc.4ysl.nc')
composite_m48('dlwr.sfc.4ysl.nc')
composite_m48('dswr.sfc.4ysl.nc')

