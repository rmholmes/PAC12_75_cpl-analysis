# /usr/bin/env python3

import numpy as np
import os
import fileinput
import glob
import sys

# Files:
experiments = ['12','14','17','20','23',
               '13','15','18','21','24',
               '26','27','28','29','30']
#'02','08','09','10','11',
               
years = ['2015','2016','2017','2018']
months = range(12)

for e in range(len(experiments)):
    for i in range(len(years)):
        for m in months:

            expstr = 'exp' + experiments[e]
            yrstr = years[i]
            mnstr = '%02d' % (m+1)

            fscr = 'fscripts/Process_' + expstr + '_' + yrstr + '_' + mnstr

            os.system('cp zonalhp_sub.sub ' + fscr)
            with fileinput.FileInput(fscr, inplace=True) as file:
                for line in file:
                    line_out = line.replace('XXEXPXX', expstr).replace('XXYRXX', yrstr).replace('XXMNXX', mnstr)
                    print(line_out, end='')
            os.system('qsub ' + fscr)
            
                    
