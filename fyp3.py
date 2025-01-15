#!/usr/bin/env python3
import subprocess
import pandas as pd
import re
import os
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM

# Read OBSID,index,z,kT and nH  from the Excel file
df = pd.read_excel('/media/sf_Portal/Clusterlist_2.xlsx')
OB=df['OBSIDs'].dropna()
N=df['Index'].dropna()
redshift=df['z'].dropna()
kT=df['kT'].dropna()
nH=df['N_H'].dropna()

# Function to calculate the norm of the cluster
def Norma(z, kT):
    H0 = 70 * u.km / u.s / u.Mpc
    Om0 = 0.3
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
    T = kT * u.keV 
    H_z = cosmo.H(z).to(u.km / u.s / u.Mpc)
    E = H_z / H0
    norm = (T**0.5) * (E**3) / ((1 + z)**4)
    return norm.value

for i,n,z,T,H in zip(OB,N,redshift,kT,nH):
    obsid_list = [item.strip() for item in str(i).split(',')]
    norm = Norma(z,T)
    command1=f'''
    source ~/ciao-4.16/bin/ciao.sh
    cd {n}/{i}
    mkdir mantz
    dmcopy "images/broad_thresh.img[exclude sky=region(sources_mod.reg)]" mantz/IM.fits clobber=yes
    dmcopy "images/broad_thresh.expmap[exclude sky=region(sources_mod.reg)]" mantz/EX.fits clobber=yes
    dmcopy "rbkg.fits[exclude sky=region(sources_mod.reg)]" mantz/BG.fits clobber=yes
    echo {z} > mantz/redshift.txt
    echo {T} > mantz/kT.txt
    echo {H} > mantz/nH.txt
    echo {norm} > mantz/norm.txt
    '''
    command2=f'''
    source ~/ciao-4.16/bin/ciao.sh
    cd {n}
    mkdir mantz
    dmcopy "sum_broad_thresh.img[exclude sky=region(sources_mod2.reg)]" mantz/IM.fits clobber=yes
    dmcopy "sum_broad_thresh.expmap[exclude sky=region(sources_mod2.reg)]" mantz/EX.fits clobber=yes
    dmcopy "sumbkg.fits[exclude sky=region(sources_mod2.reg)]" mantz/BG.fits clobber=yes
    echo {z} > mantz/redshift.txt
    echo {T} > mantz/kT.txt
    echo {H} > mantz/nH.txt
    echo {norm} > mantz/norm.txt
    '''
    if len(obsid_list)==1:
        result = subprocess.run(command1, shell=True, capture_output=True, text=True, executable='/bin/bash')
    if len(obsid_list)>1:
        result = subprocess.run(command2, shell=True, capture_output=True, text=True, executable='/bin/bash')
    print(result.stderr)
    print(f'finish {n}')   