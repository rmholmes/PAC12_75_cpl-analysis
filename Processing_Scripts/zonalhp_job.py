# /usr/bin/env python3

import numpy as np
import sys
import glob
import pac12_tools as ptools

if __name__ == "__main__":

    filt_width = 6.

    base = sys.argv[1]
    year = sys.argv[2]
    month = sys.argv[3]

    file_in_day = glob.glob(base + year + month + '*/croco_out_day.nc')[0]
    file_in_mon = file_in_day[:-6] + 'mon.nc'
    file_out = file_in_day[:-6] + 'mon_hp.nc'

    sys.stdout.write(file_in_day + ' to ' + file_out + '\n');sys.stdout.flush();
    ptools.calc_zhp_std_variables(file_in_day,file_in_mon,file_out,filt_width)

