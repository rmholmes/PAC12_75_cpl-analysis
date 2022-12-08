# /usr/bin/env python3

import numpy as np
import sys
import glob
import pac12_tools as ptools

filt_width = 12.

base = sys.argv[1]
year = sys.argv[2]

files_in = glob.glob(base + year + '*/croco_out_day.nc')
files_out = [f[:-6] + 'mon_hp.nc' for f in files_in]

for i,f_in in enumerate(files_in):
    f_out = files_out[i]
    sys.stdout.write(f_in + ' to ' + f_out + '\n')
    ptools.calc_zhp_std_variables(f_in,f_out,filt_width)

