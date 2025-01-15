#!/usr/bin/env python3
import subprocess
import pandas as pd
import re
import os
import math

RATIO=[]
RATIOERR=[]
df = pd.read_excel('/media/sf_Portal/Clusterlist_2.xlsx')
OB=df['OBSIDs'].dropna()
N=df['Index'].dropna()

for i,n in zip(OB,N):
    obsid_list = [item.strip() for item in str(i).split(',')]
    command1=f'''
    source ~/ciao-4.16/bin/ciao.sh
    cd {n}/{i}
    dmcopy "c.fits[energy=10000:12000]" h.fits clobber=yes
    dmcopy "h.fits[exclude sky=region(sources_mod.reg)]" eh.fits clobber=yes
    dmcopy "blank.evt[energy=10000:12000]" hbkg.fits clobber=yes
    dmcopy "hbkg.fits[exclude sky=region(sources_mod.reg)]" ehbkg.fits clobber=yes
    dmlist eh.fits counts
    dmlist ehbkg.fits counts
    '''

    if len(obsid_list)==1:
        result = subprocess.run(command1, shell=True, capture_output=True, text=True, executable='/bin/bash')
        matches = re.findall(r'\b\d+\b', result.stdout)
        src_counts = int(matches[-2])
        bk_counts = int(matches[-1])
        r=src_counts/bk_counts
        delta_src = math.sqrt(src_counts)
        delta_bk = math.sqrt(bk_counts)
        delta_r = math.sqrt((delta_src / bk_counts)**2 +(src_counts * delta_bk / (bk_counts**2))**2)
        command2=f'''
        source ~/ciao-4.16/bin/ciao.sh
        cd {n}/{i}
        echo -e "{r}\t{delta_r}" > mantz/bg_scaling.txt
        '''
        write = subprocess.run(command2, shell=True, capture_output=True, text=True, executable='/bin/bash')
        RATIO.append(r)
        RATIOERR.append(delta_r)

    if len(obsid_list)>1:
        command3=f'''
        source ~/ciao-4.16/bin/ciao.sh
        cd {n}
        ds9 mantz/IM.fits -smooth radius 2 -cmap b -scale log
        ds9 mantz/BG.fits -smooth radius 2 -cmap b -scale log
        '''
        rebkg = subprocess.run(command3, shell=True, capture_output=True, text=True, executable='/bin/bash')
        src_counts = float(input(print("Enter the source counts: ")))
        bkg_counts = float(input(print("Enter the background counts: ")))
        r=src_counts/bkg_counts
        delta_src = math.sqrt(src_counts)
        delta_bk = math.sqrt(bkg_counts)
        delta_r = math.sqrt((delta_src / bkg_counts)**2 +(src_counts * delta_bk / (bkg_counts**2))**2)

        command4=f'''
        source ~/ciao-4.16/bin/ciao.sh
        cd {n}
        echo -e "{r}\t{delta_r}" > mantz/bg_scaling.txt
        
        '''
        write = subprocess.run(command4, shell=True, capture_output=True, text=True, executable='/bin/bash')
        RATIO.append(r)
        RATIOERR.append(delta_r)

    print(result.stderr)
    print(f'finish {n}')
    print(RATIO)
    print(RATIOERR)