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

def zhp_filt(var,filt_width):
    # Calculate zonal filtered version of a variable

    for dim in var.dims:
        if dim.startswith('x'):
            x = dim                                # Name of zonal dimension
        else:
            raise RuntimeError("Error in zhp_filt: No zonal dimension found")

    dx = (var[x][1]-var[x][0]).values
    
    return(var.rolling({x:int(filt_width/dx)},center=True).mean())

def calc_zhp_std_variables(data_in,file_out):
    """
    Calculates:
    1) zonal-high pass filtered variances of SST, SSH and velocities
    2) Eddy and mean wind work

    from a month of daily CROCO output data and saves them back into a
    netcdf file in the same folder.
    """

    data = xr.open_dataset(data_in,chunks={'time_counter':1})

    data = create_coords_CROCO(data)

    

# Old stuff
# ---------------------------

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
