import xarray as xr
from xgcm import Grid
import numpy as np
import sys, os, glob

# General tools:
# ---------------------------------------------------------

def create_coords_CROCO(ds):

    try:
        ds["x_rho"] = ds.nav_lon_rho.isel(y_rho=0)
        ds["y_rho"] = ds.nav_lat_rho.isel(x_rho=0)
        ds["x_w"]   = ds.nav_lon_rho.isel(y_rho=0).rename({'x_rho':'x_w'})
        ds["y_w"]   = ds.nav_lat_rho.isel(x_rho=0).rename({'y_rho':'y_w'})
        ds["x_u"]   = ds.nav_lon_u.isel(y_u=0)
        ds["y_u"]   = ds.nav_lat_u.isel(x_u=0)
        ds["x_v"]   = ds.nav_lon_v.isel(y_v=0)
        ds["y_v"]   = ds.nav_lat_v.isel(x_v=0)
        ds = ds.set_coords({'x_rho','y_rho','x_u','y_u','x_v','y_v','x_w','y_w'})
    except:
        ds["x_rho"] = ds.nav_lon.isel(y=0)
        ds["y_rho"] = ds.nav_lat.isel(x=0)
        ds = ds.set_coords({'x_rho','y_rho'})

    try:
        ds["z_rho"] =  ds.z_rho.mean('time').fillna(0.)
        ds = ds.set_coords({'z_rho'})
    except:
        pass

    return(ds)

# Zonal filtering tools:
# ---------------------------------------------------------

def zlp_filt(var,filt_width):
    # Calculate zonal filtered version of a variable

    dims = var.dims
    inds = [index for index,item in enumerate(dims) if item.startswith('x')]
    if (len(inds) != 1):
        raise RuntimeError("Error in zhp_filt: less than or greater than 1 zonal dimension found")
    else:
        x = dims[inds[0]]

    dx = (var[x][1]-var[x][0]).values

    return(var.rolling({x:int(filt_width/dx)},center=True).mean())

def calc_zhp_std_variables(data_in,file_out,filt_width):
    """
    Calculates:
    1) zonal-high pass filtered variances of SST, SSH and velocities
    2) Eddy and mean wind work

    from a month of daily CROCO output data and saves them back into a
    netcdf file in the same folder.

    data_in = netcdf file of daily CROCO output data  (croco_out_day.nc)
    file_out = filename to use for output file.
    filt_width = filter window width in degrees.
    """

    data = xr.open_dataset(data_in,chunks={'time_counter':1}).rename({'time_counter':'time'})
    data = create_coords_CROCO(data)

    DT = xr.DataArray(data=len(data.time)).assign_attrs({'Name':'Number of days in averaging period'})

    # SST:
    SST_hp = data.temp_surf - zlp_filt(data.temp_surf,filt_width)
    SST_hp_var = (SST_hp**2.).mean('time').load()

    # SSH:
    SSH_hp = data.zeta - zlp_filt(data.zeta,filt_width)
    SSH_hp_var = (SSH_hp**2.).mean('time').load()

    # U, EWWU, MWWU:
    U_lp    = zlp_filt(data.u_surf,filt_width)
    U_hp    = data.u_surf - U_lp
    TAUX_lp = zlp_filt(data.sustr,filt_width)
    TAUX_hp = data.sustr - TAUX_lp

    U_hp_var = (U_hp**2.).mean('time').load()
    U_lp_var = (U_lp**2.).mean('time').load()

    EWWU = (U_hp*TAUX_hp).mean('time').load()
    MWWU = (U_lp*TAUX_lp).mean('time').load()

    # V, EWWV, MWWV:
    V_lp    = zlp_filt(data.v_surf,filt_width)
    V_hp    = data.v_surf - V_lp
    TAUY_lp = zlp_filt(data.svstr,filt_width)
    TAUY_hp = data.svstr - TAUY_lp

    V_hp_var = (V_hp**2.).mean('time').load()
    V_lp_var = (V_lp**2.).mean('time').load() # Not sure about this one...

    EWWV = (V_hp*TAUY_hp).mean('time').load()
    MWWV = (V_lp*TAUY_lp).mean('time').load()

    # Add metadata:
    SST_hp_var = SST_hp_var.assign_attrs({'Name':'SST high-pass variance'})
    SSH_hp_var = SSH_hp_var.assign_attrs({'Name':'SSH high-pass variance'})
    U_hp_var = U_hp_var.assign_attrs({'Name':'U high-pass variance'})
    V_hp_var = V_hp_var.assign_attrs({'Name':'V high-pass variance'})
    U_lp_var = U_lp_var.assign_attrs({'Name':'U low-pass variance'})
    V_lp_var = V_lp_var.assign_attrs({'Name':'V low-pass variance'})
    EWWU = EWWU.assign_attrs({'Name':'U eddy wind-work'})
    MWWU = MWWU.assign_attrs({'Name':'U mean wind-work'})
    EWWV = EWWV.assign_attrs({'Name':'V eddy wind-work'})
    MWWV = MWWV.assign_attrs({'Name':'V mean wind-work'})

    # Combine into a single Dataset and write out:
    ds = xr.Dataset(data_vars={'SSH_hp_var':SSH_hp_var,'SST_hp_var':SST_hp_var,
                  'U_hp_var':U_hp_var,'V_hp_var':V_hp_var,'U_lp_var':U_lp_var,'V_lp_var':V_lp_var,
                  'EWWU':EWWU,'MWWU':MWWU,'EWWV':EWWV,'MWWV':MWWV,'DT':DT})
    ds.to_netcdf(file_out)
                               
# Old stuff
# ----------------------------------------------------------------------------------------------------------

# Define a function to calculate depths (python script adapted from croco_tools/Preprocessing_tools/zlevs.m):
def calc_z(typ,ds):
    # Function to calculate depth of levels from CROCO output. This script has been adapted
    # from the croco_tools/Preprocessing_tools/zlevs.m script.
    # 
    # Inputs:
    # 
    # ds = CROCO history or average file xarray dataset
    # typ = 'r'; rho-points, 'w'; w-points
    #
    # Outputs:
    # z = depth of rho or w points
    
    if(ds.Vtransform != 2):
        print('ERROR - wrong Vtransform')
        return
    
    if (typ=='r'):
        Cs = ds.Cs_r
        sc = ds.sc_r
        z = xr.zeros_like(ds.temp).rename('z_rho')
        N = len(ds.s_rho)
    elif (typ=='w'):
        Cs = ds.Cs_w
        sc = ds.sc_w
        z = xr.zeros_like(ds.w).rename('z_w')
        N = len(ds.s_w)
    
    h = ds.h#.where(ds.h==0.,1.e-2,ds.h)
    #Dcrit = 0.01
    zeta = ds.zeta#.where(ds.zeta<(Dcrit-h),Dcrit-h,ds.zeta)
    
    hinv=1/h;
    h2=(h+ds.hc)
    cff = ds.hc*sc
    h2inv = 1/h2
    
    z = (cff+Cs*h)*h/h2 + zeta*(1+(cff+Cs*h)*h2inv)
    
    return(z)

# Define a function to calculate depths (python script adapted from croco_tools/Preprocessing_tools/zlevs.m):
def calc_z1D(typ,ds,h,zeta):
    # Function to calculate 1D depth of levels from CROCO output given an input H and zeta
    # 
    # Inputs:
    # 
    # ds = CROCO history or average file xarray dataset
    # typ = 'r'; rho-points, 'w'; w-points
    #
    # Outputs:
    # z = depth of rho or w points
    
    if(ds.Vtransform != 2):
        print('ERROR - wrong Vtransform')
        return
    
    if (typ=='r'):
        Cs = ds.Cs_r
        sc = ds.sc_r
        z = xr.zeros_like(Cs).rename('z_rho')
        N = len(ds.s_rho)
    elif (typ=='w'):
        Cs = ds.Cs_w
        sc = ds.sc_w
        z = xr.zeros_like(Cs).rename('z_w')
        N = len(ds.s_w)
    
    hinv=1/h;
    h2=(h+ds.hc)
    cff = ds.hc*sc
    h2inv = 1/h2
    
    z = (cff+Cs*h)*h/h2 + zeta*(1+(cff+Cs*h)*h2inv)
    
    return(z)
