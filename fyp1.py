#!/usr/bin/env python3
import subprocess
import pandas as pd
import re
import os

# Read OBSIDs and index from the Excel file
df = pd.read_excel("/media/sf_Portal/Clusterlist_2.xlsx")
OB = df["OBSIDs"].dropna()
N = df["Index"].dropna()


# Function to map ACIS-I/ACIS-S to CIAO identifiers
def replace_acis(data):
    if data == "ACIS-I":
        return "0:3"
    elif data == "ACIS-S":
        return "7"
    return data


# Iterate over OBSIDs and order values
for i, n in zip(OB, N):
    obsid_list = [item.strip() for item in str(i).split(",")]
    com = f"""
    mkdir {n}
    """
    result = subprocess.run(
        com,
        shell=True,
        capture_output=True,
        text=True,
        executable="/bin/bash",
    )
    for j in obsid_list:
        # Step 1: Find the Chandra obsid with the
        # maximum exposure time
        command_find_obsid = f"""
        source ~/ciao-4.16/bin/ciao.sh
        find_chandra_obsid {j}
        """
        result = subprocess.run(
            command_find_obsid,
            shell=True,
            capture_output=True,
            text=True,
            executable="/bin/bash",
        )
        matches = re.findall(
            r"^\s*(\d+)\s+\S+\s+(ACIS-[IS])\s+\S+\s+(\d+\.\d+)",
            result.stdout,
            re.MULTILINE,
        )
        for match in matches:
            found_obsid = int(match[0])
            instrument = match[1]
            ccdid = replace_acis(instrument)
            exposure_time = match[2]
            reduction_commands = f"""
            source ~/ciao-4.16/bin/ciao.sh
            cd {n}
            download_chandra_obsid {found_obsid}
            chandra_repro {found_obsid} outdir=''
            cd {found_obsid}
            dmcopy "repro/acisf{found_obsid :05d}_repro_evt2.fits[ccd_id={ccdid}]" c.fits clobber=yes
            fluximage c.fits images/ binsize=2 bands=broad units=area psfecf=0.9 clobber=yes
            blanksky c.fits outfile=blank.evt tmpdir=. mode=h clobber=yes
            dmcopy "blank.evt[energy=500:7000]" ebkg.fits clobber=yes
            reproject_image ebkg.fits ~/Downloads/{n}/{obsid_list[0]}/images/broad_thresh.img rbkg.fits clobber=yes
            punlearn wavdetect
            pset wavdetect regfile=sources.reg
            pset wavdetect ellsigma=4
            wavdetect images/broad_thresh.img sources.fits sources.scell sources.image sources.nbkg scales='2.0 4.0' psffile=~/Downloads/{n}/{found_obsid}/images/broad_thresh.psfmap clobber=yes
            """

            if found_obsid == int(j):
                print(
                    f"Found matching OBSID: {found_obsid}, Instrument: {instrument}, Exposure Time: {exposure_time}"
                )
                result = subprocess.run(
                    reduction_commands,
                    shell=True,
                    capture_output=True,
                    text=True,
                    executable="/bin/bash",
                )
                print(result.stderr)
            else:
                print(f"Skipping OBSID: {found_obsid}, it does not match target.")
            print(f"finish {j}")

    # For multiple OBSIDs
    if len(obsid_list) > 1:
        obsid_str = ",".join([f"{k}/rebkg.fits" for k in obsid_list])
        obsid_str2 = ",".join([f"{k}/c.fits" for k in obsid_list])
        length_str = "+".join(f"img{n}" for n in range(1, len(obsid_list) + 1))
        MERGE = f"""
        source ~/ciao-4.16/bin/ciao.sh
        cd {n}
        merge_obs {obsid_str2} sum binsize=2 bands=broad units=area psfecf=0.9 clobber=yes
        """
        mergeresult = subprocess.run(
            MERGE, shell=True, capture_output=True, text=True, executable="/bin/bash"
        )
        print(mergeresult.stderr)

        for obs in obsid_list:
            comforreprobkg = f"""
            source ~/ciao-4.16/bin/ciao.sh
            cd {n}/{obs}
            reproject_image ebkg.fits ~/Downloads/{n}/sum_broad_thresh.img rebkg.fits clobber=yes
            """
            bkgresult = subprocess.run(
                comforreprobkg,
                shell=True,
                capture_output=True,
                text=True,
                executable="/bin/bash",
            )
            print(result.stderr)
        command2 = f"""
        source ~/ciao-4.16/bin/ciao.sh
        cd {n}
        dmimgcalc {obsid_str} none sumbkg.fits op="imgout={length_str}" clobber=yes
        punlearn wavdetect
        pset wavdetect psffile=sum_broad_thresh.psfmap
        pset wavdetect regfile=sources.reg
        pset wavdetect ellsigma=4
        wavdetect sum_broad_thresh.img sources.fits sources.scell sources.image sources.nbkg scales='2.0 4.0' psffile=sum_broad_thresh.psfmap clobber=yes
        """
        result = subprocess.run(
            command2, shell=True, capture_output=True, text=True, executable="/bin/bash"
        )
        print(result.stderr)
    print(f"finish {n}")
