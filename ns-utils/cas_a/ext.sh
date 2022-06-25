#!/bin/bash -x

# This prepares HRC-S data for timing analysis

wrk_dir="/Users/silver/dat/cxo/cas/"
ra="350.86642917"
de="58.81180833"
ids=(10227 10228 10229 10892)

################################################################
# Notes
# Set the correct bad pixel list, chandra_repro creates a bad pixel
# file and sets this as the current file. So, if extraction is
# performed right after repro all's fine. However, if going back, then
# SET THE BAD PIXEL FILE TO THE ONE CREATED BY chandra_repro FOR THE
# SPECIFIC OBSERVATION
#acis_set_ardlib "acisf${ID}_repro_bpix1.fits"


################################################################

source /usr/local/ciao-4.12/bin/ciao.bash
nid=${#ids[@]}

for id in $(seq 0 $(($nid-1)))
do
    cd ${wrk_dir}${ids[id]}
    ids[id]=$(printf %05d ${ids[id]})

    # chandra_repro indir=. outdir=./repro

    cd repro
    pset ardlib AXAF_HRC-S_BADPIX_FILE=hrcf${ids[id]}_repro_bpix1.fits
    # axbary infile=hrcf${ids[id]}_repro_evt2.fits \
    #        orbitfile=$(ls ../primary/orbit*) \
    #        outfile=${ids[id]}_bc_evt.fits \
    #        ra=${ra} \
    #        dec=${de}
    dmcopy "${ids[id]}_bc_evt.fits[sky=region(src.reg)]" ${ids[id]}_bcs_evt.fits
done

echo "EOF"
