#!/usr/bin/env python3
import subprocess
import pandas as pd
import re
import os

# Read OBSIDs, index and redshift from the Excel file
df = pd.read_excel('/media/sf_Portal/Clusterlist_2.xlsx')
OB=df['OBSIDs'].dropna()
N=df['Index'].dropna()
redshift=df['z'].dropna()
S=[]
P=[]
A=[]
DS=[]
DP=[]
DA=[]
for i,n,z in zip(OB,N,redshift):
    obsid_list = [item.strip() for item in str(i).split(',')]
    print("Conducting spa calculation for", n)
    command1=f'''
    cd morph
    ./morphology.bash ~/Downloads/{n}/{i}/mantz
    ./morphology.bash ~/Downloads/{n}/{i}/mantz 1000
    '''
    command2=f'''
    cd morph
    ./morphology.bash ~/Downloads/{n}/mantz
    ./morphology.bash ~/Downloads/{n}/mantz 1000
    '''
    command3=f'''
    cd morph
    ./reduce2.R ~/Downloads/{n}/{i}/mantz/morph.log {z} ~/Downloads/{n}/{i}/mantz/morph_boot.dat
    '''
    command4=f'''
    cd morph
    ./reduce2.R ~/Downloads/{n}/mantz/morph.log {z} ~/Downloads/{n}/mantz/morph_boot.dat
    '''
    if len(obsid_list)==1:
        result1 = subprocess.run(command1, shell=True, capture_output=True, text=True, executable='/bin/bash')
        result3= subprocess.run(command3, shell=True, capture_output=True, text=True, executable='/bin/bash')
        lines = result3.stdout.strip().split('\n')
        try:
            spa_line = lines[-4]
            s= float(spa_line.split()[0])
            p= float(spa_line.split()[1])
            a= float(spa_line.split()[2])
            ds_line = lines[-3]
            dp_line = lines[-2]
            da_line = lines[-1]
            da = float(da_line.split()[2])
            dp = float(dp_line.split()[1])
            ds = float(ds_line.split()[0])
        except:
            s,p,a,ds,dp,da=0,0,0,0,0,0
            print(result3.stderr)
        S.append(s)
        P.append(p)
        A.append(a)
        DS.append(ds)
        DP.append(dp)
        DA.append(da)
    if len(obsid_list)>1:
        result2 = subprocess.run(command2, shell=True, capture_output=True, text=True, executable='/bin/bash')
        result4= subprocess.run(command4, shell=True, capture_output=True, text=True, executable='/bin/bash')
        lines = result4.stdout.strip().split('\n')
        try:
            spa_line = lines[-4]
            s= float(spa_line.split()[0])
            p= float(spa_line.split()[1])
            a= float(spa_line.split()[2])
            ds_line = lines[-3]
            dp_line = lines[-2]
            da_line = lines[-1]
            da = float(da_line.split()[2])
            dp = float(dp_line.split()[1])
            ds = float(ds_line.split()[0])
        except:
            s,p,a,ds,dp,da=0,0,0,0,0,0
            print(result4.stderr)
        S.append(s)
        P.append(p)
        A.append(a)
        DS.append(ds)
        DP.append(dp)
        DA.append(da)
    print(f'finish {n}')
print('S=',S)
print('P=',P)
print('A=',A)
print('DS=',DS)
print('DP=',DP)
print('DA=',DA)

#Write the data onto the Excel file
Array = []
Array.append(S)
Array.append(DS)
Array.append(P)
Array.append(DP)
Array.append(A)
Array.append(DA)
spa = pd.DataFrame({'Symmetry,s': Array[0], 'ds': Array[1], 'Peakiness,p': Array[2], 'dp': Array[3], 'Alignment,a': Array[4], 'da': Array[5]})
relax = pd.ExcelWriter('/media/sf_Portal/Data.xlsx')
spa.to_excel(relax)
relax.close()
print("!!!!!!FFIINNIISSHH!!!!!!")
