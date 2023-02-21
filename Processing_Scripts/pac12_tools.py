import xarray as xr
from xgcm import Grid
import numpy as np
import sys, os, glob

# Import Lisa's fortran codes:
import R_tools_fort as toolsF

# General tools:
# ---------------------------------------------------------

def create_coords_CROCO(ds):
    """
    Creates 1D coordinates from 2D coordinates assuming mercator in
    EPAC CROCO.
    """
    
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
        try:
            ds["x_rho"] = ds.nav_lon.isel(y=0)
            ds["y_rho"] = ds.nav_lat.isel(x=0)
            ds = ds.set_coords({'x_rho','y_rho'})
        except:
            pass
    try:
        ds["z_rho"] =  ds.z_rho.isel(time=0).fillna(0.)
        ds = ds.set_coords({'z_rho'})
    except:
        pass

    return(ds)

def create_coords_WRF(ds):
    """
    Creates 1D coordinates from 2D coordinates assuming mercator in
    EPAC WRF. Deal with longitude cyclicity
    """

    ds["nav_lon"] = ds.nav_lon.where(ds.nav_lon>0.,ds.nav_lon+360.)
    
    ds["x"] = ds.nav_lon.isel(y=0)
    ds["y"] = ds.nav_lat.isel(x=0)
    ds.set_coords({'x','y'})
    return(ds)

def create_coords_CROCO_TIW(ds,ds_grd):
    """
    Cleans up and creates 1D coordinates assuming mercator for
    TIW-region output in CROCO.

    Also defines some constants and grid parameters for convenience

    ds is a dataset from a _TIW.nc region netcdf file.
    ds_grd is a grid file from the parent grid.
    """

    # Rename dimensions to remove _TIW:
    obj = np.copy(ds.dims)
    red = {}
    for ob in obj:
        if '_TIW' in ob:
            red[ob] = ob.replace('_TIW','')
    ds = ds.rename(red)

    # Rename coordinates to remove _TIW:
    obj = np.copy(ds.coords)
    red = {}
    for ob in obj:
        if '_TIW' in ob:
            red[ob] = ob.replace('_TIW','')
    ds = ds.rename(red)

    # Do create_coords:
    ds = create_coords_CROCO(ds)

    # Define some constants:
    ds['rho0'] = 1025. # rho0
    L1, M1, S1, Mm, Lm, N = 60, 147, 39, 210, 1240, 36
    ds['L1'] = L1      # starting x_rho index
    ds['M1'] = M1      # starting y_rho index
    ds['S1'] = S1      # starting z_rho index

    # Extract grid variables from parent grid:
    ds['hc'] = ds_grd.hc
    ds['Cs_r'] = ds_grd.Cs_r.isel(s_rho=slice(S1,-1))
    ds['Cs_w'] = ds_grd.Cs_w.isel(s_w=slice(S1,-1))
    ds['pm'] = ds_grd.pm.isel(x=slice(L1,L1+Lm),y=slice(M1,M1+Mm)).rename({'x':'x_rho','y':'y_rho'})
    ds['pn'] = ds_grd.pn.isel(x=slice(L1,L1+Lm),y=slice(M1,M1+Mm)).rename({'x':'x_rho','y':'y_rho'})
    ds['h'] = ds_grd.h.isel(x=slice(L1,L1+Lm),y=slice(M1,M1+Mm)).rename({'x':'x_rho','y':'y_rho'})

    return(ds)

def create_WRF_xgcm(ds):
    """
    Creates a WRF xgcm grid object from WRF input dataset
    """
    # x and y dimensions/coords:
    ds = ds.drop(ds.keys())
    tmp = xr.Dataset({"tmp": (("y_v","x_u"),np.zeros((len(ds.y)-1,len(ds.x)-1)),),},
                    coords={"x_u":(ds.x[1:].values+ds.x[:-1].values)/2.,
                            "y_v": (ds.y[1:].values+ds.y[:-1].values)/2.},)
    ds = ds.assign(tmp).drop('tmp')

    wrf_grid = Grid(ds,coords={"x":{"center":"x","inner":"x_u"},
                               "y":{"center":"y","inner":"y_v"}},periodic=False)

    # Metrics:
    lat_to_m = 2*np.pi*6378000.0/360.0;
    dy_v = wrf_grid.diff(ds.nav_lat,'y')*lat_to_m
    dx_u = wrf_grid.diff(ds.nav_lon,'x')* \
                                lat_to_m*np.cos(np.pi/180*wrf_grid.interp(ds.nav_lat,'x'))
    area_psi = wrf_grid.interp(dx_u,'y')*wrf_grid.interp(dy_v,'x')
    ds = ds.assign(dy_v=dy_v).assign(dx_u=dx_u).assign(area_psi=area_psi)

    metrics = {
            ('x',): ['dx_u'], # X distances
            ('y',): ['dy_v'], # Y distances
            ('x', 'y'): ['area_psi'] # Areas
        }
    wrf_grid = Grid(ds,coords={"x":{"center":"x","inner":"x_u"},
                               "y":{"center":"y","inner":"y_v"}},periodic=False,metrics=metrics)

    return(wrf_grid)

def TIWt(ds):
    """
    Single time transpose on numpy array for use with Fortran code
    """
    return(np.transpose(ds, (2,1,0)))
        

# Zonal filtering tools:
# ---------------------------------------------------------

def zlp_filt(var,width,typ='gau',cut=0.1):
    # Calculate zonal filtered version of a variable using either:
    #
    # typ='box': a simple boxcar moving average with width of width in
    # degrees.
    #
    # typ='gau': a Gaussian filter of standard deviation width (in
    # degrees). Keep only values greater than cut in relative
    # size. This function reproduces the SST filtering used in OASIS.

    dims = var.dims
    inds = [index for index,item in enumerate(dims) if item.startswith('x')]
    if (len(inds) != 1):
        raise RuntimeError("Error in zhp_filt: less than or greater than 1 zonal dimension found")
    else:
        x = dims[inds[0]]

    dx = (var[x][1]-var[x][0]).values

    if (typ == 'gau'): # Gaussian filter
        sigx = int(width/dx)
        nn = int(3*2*sigx+1)

        xx = np.arange(0.,nn)
        dd = np.sqrt(((xx-int(nn/2))**2)/sigx**2)
        ww = 1./(np.sqrt(np.pi)*sigx)*np.exp(-dd**2)
        keepww = np.argwhere((ww.ravel() >= max(ww.ravel())*cut))
        nnok=len(keepww)
        ww = ww.ravel()[keepww.astype('int')]
        ww = ww/ww.sum()
        
        weight = xr.DataArray(ww.ravel(), dims=['window'])

        return(var.rolling({x:nnok},center=True).construct('window').dot(weight))
        
    elif (typ == 'box'): # Boxcar filter
        return(var.rolling({x:int(width/dx)},center=True).mean())
    else:
        raise RuntimeError("Error in zhp_filt: Incorrect filter type chosen. Please choose box or gau")

def calc_zhp_std_variables(file_in_day,file_in_mon,file_out,filt_width):
    """
    Calculates:
    1) zonal-high pass filtered variances of SST, SSH, velocities and surface heat fluxes
    2) Eddy and mean wind work
    3) Surface meridional heat flux V*SST
    4) Eddy and mean APE generation metrics (Bishop et al. 2020)

    from a month of daily CROCO surface output data and saves them back into a
    netcdf file in the same folder.

    file_in_day = netcdf file of daily CROCO surface output data  (croco_out_day.nc)
    file_out = filename to use for output file.
    filt_width = filter window width in degrees (see zlp_filt function above).
    """

    data = xr.open_dataset(file_in_day,chunks={'time_counter':1})
    data = create_coords_CROCO(data)

    # SST:
    SST_lp = zlp_filt(data.temp_surf,filt_width)
    SST_hp = data.temp_surf - SST_lp
    SST_hp_var = (SST_hp**2.).mean('time_counter').load()

    # SSH:
    SSH_hp = data.zeta - zlp_filt(data.zeta,filt_width)
    SSH_hp_var = (SSH_hp**2.).mean('time_counter').load()

    # U, EWWU, MWWU:
    U_lp    = zlp_filt(data.u_surf,filt_width)
    U_hp    = data.u_surf - U_lp
    TAUX_lp = zlp_filt(data.sustr,filt_width)
    TAUX_hp = data.sustr - TAUX_lp

    U_hp_var = (U_hp**2.).mean('time_counter').load()
    U_lp_var = (U_lp**2.).mean('time_counter').load()

    EWWU = (U_hp*TAUX_hp).mean('time_counter').load()
    MWWU = (U_lp*TAUX_lp).mean('time_counter').load()

    # # V, EWWV, MWWV:
    V_lp    = zlp_filt(data.v_surf,filt_width)
    V_hp    = data.v_surf - V_lp
    TAUY_lp = zlp_filt(data.svstr,filt_width)
    TAUY_hp = data.svstr - TAUY_lp

    V_hp_var = (V_hp**2.).mean('time_counter').load()
    V_lp_var = (V_lp**2.).mean('time_counter').load() # Not sure about this one...

    EWWV = (V_hp*TAUY_hp).mean('time_counter').load()
    MWWV = (V_lp*TAUY_lp).mean('time_counter').load()

    # V_hp*SST_hp for TIW meridional heat flux:
    grid = Grid(data,coords={"y":{"center":"y_rho","inner":"y_v"},"x":{"center":"x_rho","inner":"x_u"}},periodic=False)
    VSST_hp = (grid.interp(V_hp,'y').rename({'x_v':'x_rho'})*SST_hp).mean('time_counter').load()
#    USST_hp = (grid.interp(U_hp,'x').rename({'y_u':'y_rho'})*SST_hp).mean('time_counter').load() # This is not working because of some chunking error?

    # surface heat flux for APE generation:
    Q_lp = zlp_filt(data.shflx,filt_width)
    Q_hp = data.shflx - Q_lp

    Q_hp_var = (Q_hp**2.).mean('time_counter').load()
    QSST_hp  = (Q_hp*SST_hp).mean('time_counter').load()
    QSST_lp  = (Q_lp*SST_lp).mean('time_counter').load()

    # Add metadata:
    SST_hp_var = SST_hp_var.assign_attrs({'long_name':'SST high-pass variance','units':'degC^2'})
    SSH_hp_var = SSH_hp_var.assign_attrs({'long_name':'SSH high-pass variance','units':'m2'})
    U_hp_var = U_hp_var.assign_attrs({'long_name':'U high-pass variance','units':'m2s-2'})
    V_hp_var = V_hp_var.assign_attrs({'long_name':'V high-pass variance','units':'m2s-2'})
    U_lp_var = U_lp_var.assign_attrs({'long_name':'U low-pass variance','units':'m2s-2'})
    V_lp_var = V_lp_var.assign_attrs({'long_name':'V low-pass variance','units':'m2s-2'})
    EWWU = EWWU.assign_attrs({'long_name':'U eddy wind-work','units':'Nm-1s-1'})
    MWWU = MWWU.assign_attrs({'long_name':'U mean wind-work','units':'Nm-1s-1'})
    EWWV = EWWV.assign_attrs({'long_name':'V eddy wind-work','units':'Nm-1s-1'})
    MWWV = MWWV.assign_attrs({'long_name':'V mean wind-work','units':'Nm-1s-1'})
    VSST_hp = VSST_hp.assign_attrs({'long_name':'V high-pass * SST high-pass mean','units':'degCms-1'})
    Q_hp_var = Q_hp_var.assign_attrs({'long_name':'Q high-pass variance','units':'W2m-4'})
    QSST_hp = QSST_hp.assign_attrs({'long_name':'Q high-pass * SST high-pass mean','units':'WdegCm-2'})
    QSST_lp = QSST_lp.assign_attrs({'long_name':'Q low-pass * SST low-pass mean','units':'WdegC2m-2'})
#    USST = USST.assign_attrs({'long_name':'U high-pass * SST high-pass mean'})

    # Deal with time:
    DT = xr.DataArray(data=len(data.time_counter)).assign_attrs({'Name':'Number of days in averaging period'})
    data_mon = xr.open_dataset(file_in_mon,chunks={'time_counter':1})
    time_counter = data_mon.time_counter.values

    # Combine into a single Dataset and write out:
    ds = xr.Dataset(data_vars={'SSH_hp_var':SSH_hp_var,'SST_hp_var':SST_hp_var,
                  'U_hp_var':U_hp_var,'V_hp_var':V_hp_var,'U_lp_var':U_lp_var,'V_lp_var':V_lp_var,
                  'EWWU':EWWU,'MWWU':MWWU,'EWWV':EWWV,'MWWV':MWWV,'VSST_hp':VSST_hp,#'USST':USST,
                  'Q_hp_var':Q_hp_var,'QSST_hp':QSST_hp,'QSST_lp':QSST_lp,
                  'DT':DT,'time_counter':time_counter})
    
    # Add time_counter to variables:
    for v in ['SSH_hp_var','SST_hp_var','U_hp_var','U_lp_var','V_hp_var','V_lp_var','EWWU','MWWU','EWWV','MWWV','VSST_hp','Q_hp_var','QSST_hp','QSST_lp','DT']: #'USST'
        ds[v] = ds[v].expand_dims(dim={'time_counter':ds.time_counter})

    ds.encoding = {'unlimited_dims': ['time_counter']}
    ds.to_netcdf(file_out)

def calc_zhp_3d_variables(file_in_dayTIW,file_in_day,file_in_mon,file_in_grd,file_out,filt_width):
    """
    Calculates eddy correlation variables from 3D daily output in the
    TIW region for energy and heat budget quantification.

    Calculations are performed for one month of output data and then
    saved back into a netcdf file in the same folder.

    file_in_day = netcdf file of daily CROCO 3D output data in TIW region (croco_out_day_TIW.nc)
    file_out = filename to use for output file.
    filt_width = filter window width in degrees (see zlp_filt function above).
    """

    data_day = xr.open_dataset(file_in_day,chunks={'time_counter':1})
    data_dayTIW = xr.open_dataset(file_in_dayTIW,chunks={'time_counter':1})
    data_mon = xr.open_dataset(file_in_mon,chunks={'time_counter':1})
    data_grd = xr.open_dataset(file_in_grd)
    data_dayTIW = create_coords_CROCO_TIW(data_dayTIW,data_grd)

    # Setup storage arrays to store mean across months:
    arr = xr.zeros_like(data_dayTIW.temp.isel(time_counter=0)).values
    uuUx_np = np.copy(arr)
    uvUy_np = np.copy(arr)
    uvVx_np = np.copy(arr)
    vvVy_np = np.copy(arr)
    uwUz_np = np.copy(arr)
    vwVz_np = np.copy(arr)
    wb_np = np.copy(arr)

    # Deal with some coordinates:
    [Mm,Lm] = np.shape(data_dayTIW.temp.isel(time_counter=0).isel(s_rho=0))
    L1, M1, S1 = data_dayTIW['L1'].values, data_dayTIW['M1'].values, data_dayTIW['S1'].values,
    rho0 = data_dayTIW['rho0'].values

    hc = data_grd.hc
    Cs_r = data_grd.Cs_r
    Cs_w = data_grd.Cs_w
    pm = data_dayTIW.pm
    pn = data_dayTIW.pn
    h = data_dayTIW.h

    z_r = xr.zeros_like(data_mon.z_rho.isel(x_rho=slice(L1,L1+Lm),y_rho=slice(M1,M1+Mm)).isel(s_rho=slice(S1,None)).isel(time_counter=0))
    z_w = xr.zeros_like(data_dayTIW.omega.isel(time_counter=0))
    b = xr.zeros_like(data_dayTIW.temp.isel(time_counter=0))
    
    # Grid object:
    grid = Grid(data_dayTIW,coords={"x":{"center":"x_rho","outer":"x_u"},
                                    "y":{"center":"y_rho","outer":"y_v"},
                                    "s":{"center":"s_rho","outer":"s_w"}},periodic=False) #NOTE: This is different to parent grid becuase _u and _v are outer not inner.

    # start time loop:
    tL = len(data_dayTIW.time_counter.values)
    for ti in range(tL):
        print('calc_zhp_3d_variables: Doing time %03d of %03d, file ' % (ti+1,tL) + file_out)

        # Vertical grid:
        zeta = data_day.zeta.isel(x_rho=slice(L1,L1+Lm),y_rho=slice(M1,M1+Mm)).isel(time_counter=ti)
        z_rt, z_wt = toolsF.zlevs(h.T,zeta.T,hc,Cs_r,Cs_w)
        z_r.values = TIWt(z_rt[:,:,S1:])
        z_w.values = TIWt(z_wt[:,:,S1:])

        # Buoyancy calculation:
        temp = data_dayTIW.temp.isel(time_counter=ti)
        salt = data_dayTIW.salt.isel(time_counter=ti)
        b.values = TIWt(toolsF.get_buoy(TIWt(temp.values),TIWt(salt.values),TIWt(z_r.values),TIWt(z_w.values),rho0))

        u = data_dayTIW.u.isel(time_counter=ti)
        v = data_dayTIW.v.isel(time_counter=ti)
        omega = data_dayTIW.omega.isel(time_counter=ti)

        # Filtering:
        u_lp = zlp_filt(u,filt_width)
        u_hp = u - u_lp
        u_hp = u_hp.chunk({'x_u': u_hp.sizes['x_u']}) # To fix some weird apply_ufunc error below
        u_lp = u_lp.chunk({'x_u': u_lp.sizes['x_u']}) # To fix some weird apply_ufunc error below

        v_lp = zlp_filt(v,filt_width)
        v_hp = v - v_lp
        v_hp = v_hp.chunk({'x_v': v_hp.sizes['x_v']}) # To fix some weird apply_ufunc error below
        v_lp = v_lp.chunk({'x_v': v_lp.sizes['x_v']}) # To fix some weird apply_ufunc error below

        omega_lp = zlp_filt(omega,filt_width)
        omega_hp = omega - omega_lp
        omega_hp = omega_hp.chunk({'x_w': omega_hp.sizes['x_w']}) # To fix some weird apply_ufunc error below

        b_lp = zlp_filt(b,filt_width)
        b_hp = b - b_lp
        b_hp = b_hp.chunk({'x_rho': b_hp.sizes['x_rho']}) # To fix some weird apply_ufunc error below

        # Do the calculations:
        uu = grid.interp(u_hp*u_hp,'x').rename({'y_u':'y_rho'})
        vv = grid.interp(v_hp*v_hp,'y').rename({'x_v':'x_rho'})
        uv = grid.interp(u_hp,'x').rename({'y_u':'y_rho'})*grid.interp(v_hp,'y').rename({'x_v':'x_rho'})
        uw = grid.interp(u_hp,'x').rename({'y_u':'y_rho'})*grid.interp(omega_hp,'s').rename({'x_w':'x_rho','y_w':'y_rho'})
        vw = grid.interp(v_hp,'y').rename({'x_v':'x_rho'})*grid.interp(omega_hp,'s').rename({'x_w':'x_rho','y_w':'y_rho'})
        wb = b_hp*grid.interp(omega_hp,'s').rename({'x_w':'x_rho','y_w':'y_rho'})

        Ux = grid.diff(u_lp,'x').rename({'y_u':'y_rho'})*pm
        Vy = grid.diff(v_lp,'y').rename({'x_v':'x_rho'})*pn
        Uy = grid.interp(grid.interp(grid.diff(u_lp.rename({'y_u':'y_rho'}),'y'),'x'),'y')*pn
        Vx = grid.interp(grid.interp(grid.diff(v_lp.rename({'x_v':'x_rho'}),'x'),'x'),'y')*pm
        Uz = grid.interp(grid.interp(grid.diff(u_lp,'s').rename({'y_u':'y_rho'}),'x')/grid.diff(z_r,'s'),'s')
        Vz = grid.interp(grid.interp(grid.diff(v_lp,'s').rename({'x_v':'x_rho'}),'y')/grid.diff(z_r,'s'),'s')
        
        # Save the values:
        uuUx_np += (uu*Ux).values
        uvUy_np += (uv*Uy).values
        uvVx_np += (uv*Vx).values
        vvVy_np += (vv*Vy).values
        uwUz_np += (uw*Uz).values
        vwVz_np += (vw*Vz).values
        wb_np += wb.values

    # Normalize the values
    uuUx_np = uuUx_np/tL
    uvUy_np = uvUy_np/tL
    uvVx_np = uvVx_np/tL
    vvVy_np = vvVy_np/tL
    uwUz_np = uwUz_np/tL
    vwVz_np = vwVz_np/tL
    wb_np = wb_np/tL

    # Setup meta data and xarray:
    # drop_list = ['nav_lat','nav_lon','time_centered']
    arr = xr.zeros_like(data_dayTIW.temp.isel(time_counter=0)) #.drop(drop_list)
    uuUx = arr.copy().rename('uuUx').assign_attrs({'long_name':'uu*Ux','units':'m2s-3'}); uuUx.values = uuUx_np
    uvUy = arr.copy().rename('uvUy').assign_attrs({'long_name':'uv*Uy','units':'m2s-3'}); uvUy.values = uvUy_np
    uvVx = arr.copy().rename('uvVx').assign_attrs({'long_name':'uv*Vx','units':'m2s-3'}); uvVx.values = uvVx_np
    vvVy = arr.copy().rename('vvVy').assign_attrs({'long_name':'vv*Vy','units':'m2s-3'}); vvVy.values = vvVy_np
    uwUz = arr.copy().rename('uwUz').assign_attrs({'long_name':'uw*Uz','units':'m2s-3'}); uwUz.values = uwUz_np
    vwVz = arr.copy().rename('vwVz').assign_attrs({'long_name':'vw*Vz','units':'m2s-3'}); vwVz.values = vwVz_np
    wb = arr.copy().rename('wb').assign_attrs({'long_name':'wb','units':'m2s-3'}); wb.values = wb_np
    
    # Deal with time:
    DT = xr.DataArray(data=len(data_day.time_counter)).assign_attrs({'Name':'Number of days in averaging period'})
    time_counter = data_mon.time_counter.values

    # Combine into a single Dataset and write out:
    ds = xr.Dataset(data_vars={'uuUx':uuUx,
                               'uvUy':uvUy,
                               'uvVx':uvVx,
                               'vvVy':vvVy,
                               'uwUz':uwUz,
                               'vwVz':vwVz,
                               'wb':wb,
                               'DT':DT,'time_counter':time_counter})
    
    # Add time_counter to variables:
    for v in ['uuUx','uvUy','uvVx','vvVy','uwUz','vwVz','wb','DT']:
        ds[v] = ds[v].expand_dims(dim={'time_counter':ds.time_counter})

    ds.encoding = {'unlimited_dims': ['time_counter']}
    ds.to_netcdf(file_out)

def calc_zhp_3dstd_variables(file_in_dayTIW,file_in_day,file_in_mon,file_in_grd,file_out,filt_width):
    """
    Calculates some more eddy correlation variables from 3D daily output in the
    TIW region for energy and heat budget quantification.

    Calculations are performed for one month of output data and then
    saved back into a netcdf file in the same folder.

    file_in_day = netcdf file of daily CROCO 3D output data in TIW region (croco_out_day_TIW.nc)
    file_out = filename to use for output file.
    filt_width = filter window width in degrees (see zlp_filt function above).
    """

    data_day = xr.open_dataset(file_in_day,chunks={'time_counter':1})
    data_dayTIW = xr.open_dataset(file_in_dayTIW,chunks={'time_counter':1})
    data_mon = xr.open_dataset(file_in_mon,chunks={'time_counter':1})
    data_grd = xr.open_dataset(file_in_grd)
    data_dayTIW = create_coords_CROCO_TIW(data_dayTIW,data_grd)

    # Setup storage arrays to store mean across months:
    arr = xr.zeros_like(data_dayTIW.temp.isel(time_counter=0)).values
    uu_np = np.copy(arr)
    vv_np = np.copy(arr)
    uT_np = np.copy(arr)
    vT_np = np.copy(arr)
    wT_np = np.copy(arr)

    # Grid object:
    grid = Grid(data_dayTIW,coords={"x":{"center":"x_rho","outer":"x_u"},
                                    "y":{"center":"y_rho","outer":"y_v"},
                                    "s":{"center":"s_rho","outer":"s_w"}},periodic=False) #NOTE: This is different to parent grid becuase _u and _v are outer not inner.

    # start time loop:
    tL = len(data_dayTIW.time_counter.values)
    for ti in range(tL):
        print('calc_zhp_3d_variables: Doing time %03d of %03d, file ' % (ti+1,tL) + file_out)

        u = data_dayTIW.u.isel(time_counter=ti)
        v = data_dayTIW.v.isel(time_counter=ti)
        omega = data_dayTIW.omega.isel(time_counter=ti)
        temp = data_dayTIW.temp.isel(time_counter=ti)

        # Filtering:
        u_lp = zlp_filt(u,filt_width)
        u_hp = u - u_lp
        u_hp = u_hp.chunk({'x_u': u_hp.sizes['x_u']}) # To fix some weird apply_ufunc error below
        u_lp = u_lp.chunk({'x_u': u_lp.sizes['x_u']}) # To fix some weird apply_ufunc error below

        v_lp = zlp_filt(v,filt_width)
        v_hp = v - v_lp
        v_hp = v_hp.chunk({'x_v': v_hp.sizes['x_v']}) # To fix some weird apply_ufunc error below
        v_lp = v_lp.chunk({'x_v': v_lp.sizes['x_v']}) # To fix some weird apply_ufunc error below

        omega_lp = zlp_filt(omega,filt_width)
        omega_hp = omega - omega_lp
        omega_hp = omega_hp.chunk({'x_w': omega_hp.sizes['x_w']}) # To fix some weird apply_ufunc error below

        T_lp = zlp_filt(temp,filt_width)
        T_hp = temp - T_lp
        T_hp = T_hp.chunk({'x_rho': T_hp.sizes['x_rho']}) # To fix some weird apply_ufunc error below

        # Do the calculations:
        uu = grid.interp(u_hp*u_hp,'x').rename({'y_u':'y_rho'})
        vv = grid.interp(v_hp*v_hp,'y').rename({'x_v':'x_rho'})
        uT = grid.interp(u_hp,'x').rename({'y_u':'y_rho'})*T_hp
        vT = grid.interp(v_hp,'y').rename({'x_v':'x_rho'})*T_hp
        wT = T_hp*grid.interp(omega_hp,'s').rename({'x_w':'x_rho','y_w':'y_rho'})
        
        # Save the values:
        uu_np += uu.values
        vv_np += vv.values
        uT_np += uT.values
        vT_np += vT.values
        wT_np += wT.values

    # Normalize the values
    uu_np = uu_np/tL
    vv_np = vv_np/tL
    uT_np = uT_np/tL
    vT_np = vT_np/tL
    wT_np = wT_np/tL

    # Setup meta data and xarray:
    # drop_list = ['nav_lat','nav_lon','time_centered']
    arr = xr.zeros_like(data_dayTIW.temp.isel(time_counter=0)) #.drop(drop_list)
    uu = arr.copy().rename('uu').assign_attrs({'long_name':'uu','units':'m2s-2'}); uu.values = uu_np
    vv = arr.copy().rename('vv').assign_attrs({'long_name':'vv','units':'m2s-2'}); vv.values = vv_np
    uT = arr.copy().rename('uT').assign_attrs({'long_name':'uT','units':'ms-1degC'}); uT.values = uT_np
    vT = arr.copy().rename('vT').assign_attrs({'long_name':'vT','units':'ms-1degC'}); vT.values = vT_np
    wT = arr.copy().rename('wT').assign_attrs({'long_name':'wT','units':'ms-1degC'}); wT.values = wT_np
    
    # Deal with time:
    DT = xr.DataArray(data=len(data_day.time_counter)).assign_attrs({'Name':'Number of days in averaging period'})
    time_counter = data_mon.time_counter.values

    # Combine into a single Dataset and write out:
    ds = xr.Dataset(data_vars={'uu':uu,
                               'vv':vv,
                               'uT':uT,
                               'vT':vT,
                               'wT':wT,
                               'DT':DT,'time_counter':time_counter})
    
    # Add time_counter to variables:
    for v in ['uu','vv','uT','vT','wT','DT']:
        ds[v] = ds[v].expand_dims(dim={'time_counter':ds.time_counter})

    ds.encoding = {'unlimited_dims': ['time_counter']}
    ds.to_netcdf(file_out)

def calc_zhp_wrf_variables(file_in_day,file_in_mon,file_out,filt_width):
    """

    Calculates zonal-high-pass filtered variances of WRF surface
    variables from WRF daily data and saves back into a monthly netcdf
    file in the same folder.

    file_in_day = netcdf file of daily WRF surface output data  (wrf3d_1D_*.nc)
    file_in_month = netcdf file for monthly CROCO output (just to get time value)
    file_out = filename to use for output file.
    filt_width = filter window width in degrees (see zlp_filt function above).
    """

    data = xr.open_dataset(file_in_day,chunks={'time_counter':1})
    data = create_coords_WRF(data)
    wrf_grid = create_WRF_xgcm(data)

    # SST:
    SST_hp = data.SST - zlp_filt(data.SST,filt_width)

    # Surface winds:
    WX_hp = data.U_PHYL1 - zlp_filt(data.U_PHYL1,filt_width)
    WY_hp = data.V_PHYL1 - zlp_filt(data.V_PHYL1,filt_width)

    # Surface stresses:
    SX_hp = data.TAUX - zlp_filt(data.TAUX,filt_width)
    SY_hp = data.TAUY - zlp_filt(data.TAUY,filt_width)

    # Fix chunking error:
    WX_hp = WX_hp.chunk({'x': WX_hp.sizes['x']})
    WY_hp = WY_hp.chunk({'x': WY_hp.sizes['x']})
    SX_hp = SX_hp.chunk({'x': SX_hp.sizes['x']})
    SY_hp = SY_hp.chunk({'x': SY_hp.sizes['x']})

    # Calculate divergence and curl:
    Wdiv  = wrf_grid.interp(wrf_grid.derivative(WX_hp,'x'),'y')+wrf_grid.interp(wrf_grid.derivative(WY_hp,'y'),'x')
    Wcur = wrf_grid.interp(wrf_grid.derivative(WY_hp,'x'),'y')-wrf_grid.interp(wrf_grid.derivative(WX_hp,'y'),'x')
    Sdiv  = wrf_grid.interp(wrf_grid.derivative(SX_hp,'x'),'y')+wrf_grid.interp(wrf_grid.derivative(SY_hp,'y'),'x')
    Scur = wrf_grid.interp(wrf_grid.derivative(SY_hp,'x'),'y')-wrf_grid.interp(wrf_grid.derivative(SX_hp,'y'),'x')

    # Calculate magnitudes:
    Wmag = np.sqrt((WX_hp)**2.+(WY_hp)**2.)
    Smag = np.sqrt((SX_hp)**2.+(SY_hp)**2.)

    # Do variance calculations:
    Wdiv_hp_var = (Wdiv**2.).mean('time_counter').load()
    Wcur_hp_var = (Wcur**2.).mean('time_counter').load()
    Sdiv_hp_var = (Sdiv**2.).mean('time_counter').load()
    Scur_hp_var = (Scur**2.).mean('time_counter').load()
    WX_hp_var   = (WX_hp**2.).mean('time_counter').load()
    WY_hp_var   = (WY_hp**2.).mean('time_counter').load()
    SX_hp_var   = (SX_hp**2.).mean('time_counter').load()
    SY_hp_var   = (SY_hp**2.).mean('time_counter').load()
    SST_hp_var  = (SST_hp**2.).mean('time_counter').load()
    
    # Add metadata:
    SST_hp_var  = SST_hp_var.assign_attrs({'long_name':'SST high-pass variance','units':'degC2'})
    Wdiv_hp_var = Wdiv_hp_var.assign_attrs({'long_name':'high-pass wind divergence variance','units':'s-2'})
    Wcur_hp_var = Wcur_hp_var.assign_attrs({'long_name':'high-pass wind curl variance','units':'s-2'})
    Sdiv_hp_var = Sdiv_hp_var.assign_attrs({'long_name':'high-pass wind divergence variance','units':'N2m-6'})
    Scur_hp_var = Scur_hp_var.assign_attrs({'long_name':'high-pass wind curl variance','units':'N2m-6'})
    WX_hp_var   = WX_hp_var.assign_attrs({'long_name':'high-pass U-wind variance','units':'m2s-2'})
    WY_hp_var   = WY_hp_var.assign_attrs({'long_name':'high-pass V-wind variance','units':'m2s-2'})
    SX_hp_var   = SX_hp_var.assign_attrs({'long_name':'high-pass U-stress variance','units':'N2m-4'})
    SY_hp_var   = SY_hp_var.assign_attrs({'long_name':'high-pass V-stress variance','units':'N2m-4'})

    # Deal with time:
    DT = xr.DataArray(data=len(data.time_counter)).assign_attrs({'Name':'Number of days in averaging period'})
    data_mon = xr.open_dataset(file_in_mon,chunks={'time_counter':1})
    time_counter = data_mon.time_counter.values

    # Combine into a single Dataset and write out:
    ds = xr.Dataset(data_vars={'SST_hp_var':SST_hp_var,
                               'Wdiv_hp_var':Wdiv_hp_var,
                               'Wcur_hp_var':Wcur_hp_var,
                               'Sdiv_hp_var':Sdiv_hp_var,
                               'Scur_hp_var':Scur_hp_var,
                               'WX_hp_var':WX_hp_var,
                               'WY_hp_var':WY_hp_var,
                               'SX_hp_var':SX_hp_var,
                               'SY_hp_var':SY_hp_var,
                               'DT':DT,'time_counter':time_counter})
    
    # Add time_counter to variables:
    for v in ['SST_hp_var','Wdiv_hp_var','Wcur_hp_var','Sdiv_hp_var','Scur_hp_var','WX_hp_var','WY_hp_var','SX_hp_var','SY_hp_var']:
        ds[v] = ds[v].expand_dims(dim={'time_counter':ds.time_counter})

    ds.encoding = {'unlimited_dims': ['time_counter']}
    ds.to_netcdf(file_out)

# # Old stuff
# # ----------------------------------------------------------------------------------------------------------

# # Define a function to calculate depths (python script adapted from croco_tools/Preprocessing_tools/zlevs.m):

# def calc_z(typ,ds):
#     # Function to calculate depth of levels from CROCO output. This script has been adapted
#     # from the croco_tools/Preprocessing_tools/zlevs.m script.
#     # 
#     # Inputs:
#     # 
#     # ds = CROCO history or average file xarray dataset
#     # typ = 'r'; rho-points, 'w'; w-points
#     #
#     # Outputs:
#     # z = depth of rho or w points
    
#     if(ds.Vtransform != 2):
#         print('ERROR - wrong Vtransform')
#         return
    
#     if (typ=='r'):
#         Cs = ds.Cs_r
#         sc = ds.sc_r
#         z = xr.zeros_like(ds.temp).rename('z_rho')
#         N = len(ds.s_rho)
#     elif (typ=='w'):
#         Cs = ds.Cs_w
#         sc = ds.sc_w
#         z = xr.zeros_like(ds.w).rename('z_w')
#         N = len(ds.s_w)
    
#     h = ds.h#.where(ds.h==0.,1.e-2,ds.h)
#     #Dcrit = 0.01
#     zeta = ds.zeta#.where(ds.zeta<(Dcrit-h),Dcrit-h,ds.zeta)
    
#     hinv=1/h;
#     h2=(h+ds.hc)
#     cff = ds.hc*sc
#     h2inv = 1/h2
    
#     z = (cff+Cs*h)*h/h2 + zeta*(1+(cff+Cs*h)*h2inv)
    
#     return(z)

# # Define a function to calculate depths (python script adapted from croco_tools/Preprocessing_tools/zlevs.m):
# def calc_z1D(typ,ds,h,zeta):
#     # Function to calculate 1D depth of levels from CROCO output given an input H and zeta
#     # 
#     # Inputs:
#     # 
#     # ds = CROCO history or average file xarray dataset
#     # typ = 'r'; rho-points, 'w'; w-points
#     #
#     # Outputs:
#     # z = depth of rho or w points
    
#     if(ds.Vtransform != 2):
#         print('ERROR - wrong Vtransform')
#         return
    
#     if (typ=='r'):
#         Cs = ds.Cs_r
#         sc = ds.sc_r
#         z = xr.zeros_like(Cs).rename('z_rho')
#         N = len(ds.s_rho)
#     elif (typ=='w'):
#         Cs = ds.Cs_w
#         sc = ds.sc_w
#         z = xr.zeros_like(Cs).rename('z_w')
#         N = len(ds.s_w)
    
#     hinv=1/h;
#     h2=(h+ds.hc)
#     cff = ds.hc*sc
#     h2inv = 1/h2
    
#     z = (cff+Cs*h)*h/h2 + zeta*(1+(cff+Cs*h)*h2inv)
    
#     return(z)
