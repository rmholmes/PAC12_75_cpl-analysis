# /usr/bin/env python3

import numpy as np
import os
import fileinput
import glob
import sys

# Files:
experiments = ['26']#,'09','17','18','19']
# experiments = ['02','05','06','08','09','10','11','12','13',
#               '15','16','17','18','19','20','21','24','25']
# experiments = ['22','23']
years = ['2015']#,'2016','2017','2018']
months = [3,4]#range(12)

for e in range(len(experiments)):
    for i in range(len(years)):
        for m in months:

            expstr = 'exp' + experiments[e]
            yrstr = years[i]
            mnstr = '%02d' % (m+1)

            fscr = 'fscripts/Process_' + expstr + '_' + yrstr + '_' + mnstr

            # 3d vars:
            os.system('cp zonalhp_3d_sub.sub ' + fscr)
            with fileinput.FileInput(fscr, inplace=True) as file:
                for line in file:
                    line_out = line.replace('XXEXPXX', expstr).replace('XXYRXX', yrstr).replace('XXMNXX', mnstr)
                    print(line_out, end='')
            os.system('qsub ' + fscr)
            
                    
