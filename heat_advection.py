import numpy as np
import xarray as xr
import sys
sys.path.append('//ccc/cont005/home/obsmip/maillali/TOOLS/Tracer_balance_code_LMaillard_v2/')
import R_tools_fort as toolsF


# To be looped on each time step
# Warning ! This script takes the fields u,v,temp as [lon X lat X s_rho] and zeta as [lon X lat] so the fields need to be transposed.
# Cs_r and Cs_w can normally be obtained in the croco output netcdf history. 
# pm,pn,h,hc can usually be found in the grid file 
# dt_model is the time step of the model (not sure wether it is really used here or not). 

print('... compute zlevs')
[zr,zw] = toolsF.zlevs(h,zeta, hc, Cs_r,Cs_w)
print('... compute omega')
omega = toolsF.get_omega(u,v,zr,zw,pm,pn)
print('... compute advection')
TXadv,TYadv,TVadv = toolsF.get_tracer_evolution_adv2(u,v,zr,zw,pm,pn,dt_model,t,omega)



