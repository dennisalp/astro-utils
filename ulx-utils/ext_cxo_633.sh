#!/bin/bash

function get_img {
    # Apply a 0.3 to 10 keV filter
    dmcopy "acisf${ids[id]}N003_evt2.fits[energy=300:10000]" \
	   "acisf${ids[id]}N003_evt2_03_10.fits"

    # Create a background light curve
    dmextract "acisf${ids[id]}N003_evt2_03_10.fits[bin time=::2048]" \
    	      "acisf${ids[id]}N003_bkg_lc.fits" \
    	      "opt=ltc1"

    # Find the GTI, uses clean instead of sigma
    deflare "acisf${ids[id]}N003_bkg_lc.fits" \
    	    "acisf${ids[id]}N003_gti.fits" \
    	    "method=sigma" \
	    "stddev=3" \
    	    "save=acisf${ids[id]}N003_bkg_lc.ps"

    # Apply the GTI
    dmcopy "acisf${ids[id]}N003_evt2_03_10.fits[@acisf${ids[id]}N003_gti.fits]" \
           "acisf${ids[id]}N003_evt2_03_10_clean.fits"

    dmlist acisf${ids[id]}N003_evt2_03_10.fits header | grep EXPO
    dmlist acisf${ids[id]}N003_evt2_03_10_clean.fits header | grep EXPO
    
    # ds9 print regions
    ds9 acisf${ids[id]}N003_evt2_03_10_clean.fits -scale log -cmap Heat -regions ../../src_${ids[id]}.reg -regions ../../bkg_${ids[id]}.reg -bin buffersize 2048 -zoom 5 -scale linear -pan to ${coords} wcs -print destination file -print filename acisf${ids[id]}N003_evt2_03_10_clean_img.ps -print -exit
}

# This gets one zeroth order spectrum for the entire source (i.e. standard)
# refcoord expects de_deg to include the sign, so "+" needs to be added if de_deg isn't negative
function get_spe {
    specextract infile="acisf${ids[id]}N003_evt2_03_10_clean.fits[sky=region(../../src_${ids[id]}.reg)]" \
                bkgfile="acisf${ids[id]}N003_evt2_03_10_clean.fits[sky=region(../../bkg_${ids[id]}.reg)]" \
    		outroot="acisf${ids[id]}N003_0th" \
    		correctpsf=yes \
    		weight=no \
                grouptype=NUM_CTS \
                binspec=20
    
    # Rename and link the files
    fthedit "acisf${ids[id]}N003_0th_grp.pi[0]" RESPFILE add "acisf${ids[id]}N003_0th.rmf"
    fthedit "acisf${ids[id]}N003_0th_grp.pi[0]" ANCRFILE add "acisf${ids[id]}N003_0th.corr.arf"
    fthedit "acisf${ids[id]}N003_0th_grp.pi[1]" RESPFILE add "acisf${ids[id]}N003_0th.rmf"
    fthedit "acisf${ids[id]}N003_0th_grp.pi[1]" ANCRFILE add "acisf${ids[id]}N003_0th.corr.arf"
    fthedit "acisf${ids[id]}N003_0th_grp.pi[0]" BACKFILE add "acisf${ids[id]}N003_0th_bkg.pi"
    fthedit "acisf${ids[id]}N003_0th_grp.pi[1]" BACKFILE add "acisf${ids[id]}N003_0th_bkg.pi"
}

function clean {
    for ff in *.ps
    do
	ps2pdf ${ff}
    done
    rm *.ps
}

################################################################

wrk_dir="/Users/silver/dat/cxo/633/"
ids=(2901 4536)
# ids=(4536)
coords="00:41:25.870 +40:58:45.70"



source /usr/local/ciao-4.12/bin/ciao.bash
nid=${#ids[@]}

for id in $(seq 0 $(($nid-1)))
do
    cd ${wrk_dir}${ids[id]}
    ids[id]=$(printf %05d ${ids[id]})
    echo ${ids[id]}

    # chandra_repro indir=. outdir=./repro

    cd primary
    acis_set_ardlib acisf${ids[id]}N003_bpix1.fits
    get_img
    get_spe
    clean
done

echo "EOF"
