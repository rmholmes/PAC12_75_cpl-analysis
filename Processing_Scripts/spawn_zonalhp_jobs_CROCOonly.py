# /usr/bin/env python3

import numpy as np
import os
import fileinput
import glob
import sys

# Files:
experiments = ['11']
years = ['2014','2015','2016','2017','2018']

for e in range(len(experiments)):
    for i in range(len(years)):

        expstr = 'exp' + experiments[e]
        yrstr = years[i]

        fscr = 'fscripts/ProcessCROCOonly_' + expstr + '_' + yrstr

        os.system('cp zonalhp_sub_CROCOonly.sub ' + fscr)
        with fileinput.FileInput(fscr, inplace=True) as file:
            for line in file:
                line_out = line.replace('XXEXPXX', expstr).replace('XXYRXX', yrstr)
                print(line_out, end='')
        os.system('qsub ' + fscr)
            
                    
