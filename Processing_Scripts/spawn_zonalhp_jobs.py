# /usr/bin/env python3

import numpy as np
import os
import fileinput
import glob

# Files:
experiments = ['14'] # Testing on 14.
# experiments = ['02','05','06','08','09','10','11','12','13',
#               '15','16','17','18','19','20','21,'22','23','24','25']
years = ['2015','2016','2017','2018']
# experiments = ['13']#,'06']#,'08','09','10','11','12','13']
# years = ['2018']#,'2016']#,'2017','2018']

for e in range(len(experiments)):
    for i in range(len(years)):

        expstr = 'exp' + experiments[e]
        yrstr = years[i]

        fscr = 'fscripts/Process_' + expstr + '_' + yrstr
        os.system('cp zonalhp_sub.sub ' + fscr)
        with fileinput.FileInput(fscr, inplace=True) as file:
            for line in file:
                line_out = line.replace('XXEXPXX', expstr).replace('XXYRXX', yrstr)
                print(line_out, end='')
        os.system('qsub ' + fscr)
            
                    
