#!/bin/bash
# kirsch06
# https://heasarc.gsfc.nasa.gov/docs/xmm/abc/node10.html

id="0500880201"
wrk="/Users/silver/dat/xmm/cx1/${id}_repro"
dat="/Users/silver/dat/xmm/cx1/${id}"
ra="299.59031592"
de="35.20160611"

if [[ $wrk =~ " " ]]
then
    echo "Path to working is not allowed to contain spaces! (SAS issue)"
    exit 1
fi
mkdir -p ${wrk}
cd ${wrk}

export SAS_DIR="/Users/silver/sas_18.0.0-Darwin-16.7.0-64/xmmsas_20190531_1155"
export SAS_CCFPATH="/Users/silver/ccf"
export SAS_ODF="${dat}/ODF"
export SAS_CCF="${wrk}/ccf.cif"
. ${SAS_DIR}/setsas.sh

# cifbuild
export SAS_CCF=ccf.cif
# odfingest
export SAS_ODF=$(ls -1 *SUM.SAS)

# Raw event list
# epproc burst=yes
# mv *TimingEvts.ds epn.evt

# Image
# evselect table=epn.evt \
#          withimageset=yes \
#          imageset=image.fits \
#          xcolumn=RAWX \
#          ycolumn=RAWY \
#          imagebinning=binSize \
#          ximagebinsize=1 \
#          yimagebinsize=1

# Standard filters
# evselect table=epn.evt \
#          withfilteredset=yes \
#          expression='(PATTERN <= 4)&&(PI in [200:15000])&&#XMMEA_EP' \
#          filteredset=epn_cl.evt \
#          filtertype=expression \
#          keepfilteroutput=yes \
#          updateexposure=yes \
#          filterexposure=yes

# Light curve for flares
# evselect table=epn_cl.evt \
#          withrateset=yes \
#          rateset=bkg.lc \
#          maketimecolumn=yes \
#          timecolumn=TIME \
#          timebinsize=50 \
#          makeratecolumn=yes
# dsplot table=bkg.lc x=TIME y=RATE

# GTI
# tabgtigen table=bkg.lc \
# 	  expression='RATE<=100' \
# 	  gtiset=epn.gti

# Clean event list
# evselect table=epn_cl.evt \
# 	 withfilteredset=yes \
# 	 expression="(PATTERN <= 12) && (PI in [300:10000]) && #XMMEA_EP && gti(epn.gti,TIME) && (RAWX in [27:47])" \
# 	 filteredset=epn_vcl.evt \
# 	 filtertype=expression \
# 	 keepfilteroutput=yes \
# 	 updateexposure=yes \
# 	 filterexposure=yes

# Clean image
# evselect table=epn_vcl.evt \
#          withimageset=yes \
#          imageset=image_vcl.fits \
#          xcolumn=RAWX \
#          ycolumn=RAWY \
#          imagebinning=binSize \
#          ximagebinsize=1 \
#          yimagebinsize=1

# Check for pile-up
# epatplot set=epn_vcl.evt \
#          plotfile=pn_epat.ps \
#          useplotfile=yes

# Barycentric correction
cp epn_vcl.evt epn_bc.evt
barycen table=epn_bc.evt:EVENTS \
        ephemeris=DE405 \
	srcra=${ra} \
	srcdec=${de}
