import xarray as xr
from xgcm import Grid
import numpy as np
import sys, os, glob

# General tools:
# ---------------------------------------------------------

def do_coords(ds):

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

    d["x_rho"] = d.nav_lon_rho.isel(y_rho=0)
    d["y_rho"] = d.nav_lat_rho.isel(x_rho=0)
    d["x_w"]   = d.nav_lon_rho.isel(y_rho=0).rename({'x_rho':'x_w'})
    d["y_w"]   = d.nav_lat_rho.isel(x_rho=0).rename({'y_rho':'y_w'})
    d["x_u"]   = d.nav_lon_u.isel(y_u=0)
    d["y_u"]   = d.nav_lat_u.isel(x_u=0)
    d["x_v"]   = d.nav_lon_v.isel(y_v=0)
    d["y_v"]   = d.nav_lat_v.isel(x_v=0)
    d = d.set_coords({'x_rho','y_rho','x_u','y_u','x_v','y_v','x_w','y_w'})
    


    
# Files:
PI_or_his = 1
obase = '/scratch/e14/rmh561/access-cm2/HCvar/'

if PI_or_his==1:
    base = '/g/data/p66/cm2704/archive/bi889/history/ocn/';
    name = 'PIcontrolTb05'
elif PI_or_his==0:
    base = '/g/data/p66/cm2704/archive/bj594/history/ocn/';
    name = 'hisTb05'
elif PI_or_his==2:
    base = '/g/data/p66/cm2704/archive/by350/history/ocn/';
    name = 'hisNATe1Tb05'
elif PI_or_his==3:
    base = '/g/data/p66/cm2704/archive/by438/history/ocn/';
    name = 'hisNATe2'
elif PI_or_his==4:
    base = '/g/data/p66/cm2704/archive/by563/history/ocn/';
    name = 'hisNATe3'
elif PI_or_his==5:
    base = '/g/data/p66/cm2704/archive/bw966/history/ocn/';
    name = 'hisAERe1'
elif PI_or_his==6:
    base = '/g/data/p66/cm2704/archive/bu010/history/ocn/';
    name = 'hisGHGe1'

files = glob.glob(base + 'ocean_month.nc-*')
dates = [x[-8:] for x in files]

dates = sorted(dates,key=lambda x: int(x))

# Mask:
#msk = 'ypm60'
msk = ''

# Split into sets of 600 runs for PI control (ss is 0,1):
ss = 1
dates = dates[ss*600:ss*600+600]
# dates = ['11500630']
# dates = ['19640630','18591231']

# Testing (do only 3):
# dates = dates[:3]

print(dates)
print(len(dates))#len(dates))
for i in range(len(dates)):
    fname = base + 'ocean_month.nc-' + dates[i]
    oname = obase + 'CM2_' + name + '_' + msk + '_' + dates[i] + '.mat'
    fscr = 'fscripts/Process_ACCESS_' + msk + '_' + dates[i]
    os.system('cp Process_ACCESS_CM2 ' + fscr)
    with fileinput.FileInput(fscr, inplace=True) as file:
        for line in file:
            line_out = line.replace('XXNAMEXX', 'P' + dates[i]).replace('XXFNAMEXX', fname).replace('XXONAMEXX', oname).replace('XXMSKXX', msk)
            print(line_out, end='')
    os.system('qsub ' + fscr)
            
                    
