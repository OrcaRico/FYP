#!/usr/bin/env python3
import subprocess
import pandas as pd
import re
import os

# Read OBSIDs and index from the Excel file
df = pd.read_excel('/media/sf_Portal/Clusterlist_2.xlsx')
OB=df['OBSIDs'].dropna()
N=df['Index'].dropna()
for i,n in zip(OB,N):
    obsid_list = [item.strip() for item in str(i).split(',')]
    command1=f'''
    source ~/ciao-4.16/bin/ciao.sh
    cd {n}/{i}
    ds9 images/broad_thresh.img -smooth radius 2 -cmap b -region sources.reg -scale log
    '''
    command2=f'''                                                                                                                                                         
    source ~/ciao-4.16/bin/ciao.sh                                                                                                                                       
    cd {n}
    ds9 sum_broad_thresh.img -smooth radius 2 -cmap b -region sources.reg -scale log                                                                                     
    '''
    if len(obsid_list)==1:
        result = subprocess.run(command1, shell=True, capture_output=True, text=True, executable='/bin/bash')
    if len(obsid_list)>1:
        result = subprocess.run(command2, shell=True, capture_output=True, text=True, executable='/bin/bash')
    print(result.stderr)
    print(f'finish {n}') 