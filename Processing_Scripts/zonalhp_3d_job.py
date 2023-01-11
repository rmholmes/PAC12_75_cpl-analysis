# /usr/bin/env python3

import numpy as np
import sys
import glob
import pac12_tools as ptools

filt_width = 6.

base = sys.argv[1]
year = sys.argv[2]

files_in_day = glob.glob(base + year + '*/croco_out_day.nc')
files_in_dayTIW = [f[:-6] + 'day_TIW.nc' for f in files_in_day]
files_in_mon = [f[:-6] + 'mon.nc' for f in files_in_day]
files_in_grd = [f[:-6] + 'grd.nc' for f in files_in_day]
files_out = [f[:-6] + 'mon_TIWhp.nc' for f in files_in_day]

for i,f_in in enumerate(files_in_day):
    f_in_dayTIW = files_in_dayTIW[i]
    f_in_grd = files_in_grd[i]
    f_in_mon = files_in_mon[i]
    f_out = files_out[i]
    sys.stdout.write(f_in + ' to ' + f_out + '\n')
    ptools.calc_zhp_3d_variables(f_in_dayTIW,f_in,f_in_mon,f_in_grd,f_out,filt_width)

