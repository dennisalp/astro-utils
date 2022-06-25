#!/bin/bash

id="0822580401"
wrk="/Users/silver/dat/xmm/cow/${id}_repro"
dat="/Users/silver/dat/xmm/cow/${id}"
src="circle(25969.994,28042.520,400)"
bkg="circle(28027.063,24081.554,640)"
ra="25969.994"
de="28042.520"

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
# epproc
# mv *ImagingEvts.ds epn.evt

# Light curve for flares
# evselect table=epn.evt \
# 	 withrateset=Y \
# 	 rateset=bkg.lc \
# 	 maketimecolumn=Y \
# 	 timebinsize=100 \
# 	 makeratecolumn=Y \
# 	 expression='#XMMEA_EP && (PI>10000&&PI<12000) && (PATTERN==0)'
# dsplot table=bkg.lc x=TIME y=RATE

# GTI
# tabgtigen table=bkg.lc \
# 	  expression='RATE<=0.6' \
# 	  gtiset=epn.gti

# Clean event list
# evselect table=epn.evt \
# 	 withfilteredset=yes \
# 	 expression="(PATTERN <= 12) && (PI in [150:15000]) && #XMMEA_EP && gti(epn.gti,TIME)" \
# 	 filteredset=epn_cl.evt \
# 	 filtertype=expression \
# 	 keepfilteroutput=yes \
# 	 updateexposure=yes \
# 	 filterexposure=yes

# Clean image
# evselect table=epn_cl.evt \
# 	 withimageset=yes \
# 	 imageset=epn_cl.img \
# 	 xcolumn=X \
# 	 ycolumn=Y \
# 	 imagebinning=imageSize \
#          ximagebinsize=80 \
#          yimagebinsize=80

# Source events
# evselect table=epn_cl.evt \
#          withfilteredset=yes \
# 	 expression="(X,Y) IN ${src}" \
# 	 filteredset=epn_cl_src.evt \
# 	 filtertype=expression \
# 	 keepfilteroutput=yes \
# 	 updateexposure=yes \
# 	 filterexposure=yes

# Control events
# evselect table=epn_cl.evt \
#          withfilteredset=yes \
# 	 expression="(X,Y) IN ${bkg}" \
# 	 filteredset=epn_cl_bkg.evt \
# 	 filtertype=expression \
# 	 keepfilteroutput=yes \
# 	 updateexposure=yes \
# 	 filterexposure=yes

# Clean light curve
# evselect table=epn_cl_src.evt \
# 	 withrateset=yes \
# 	 rateset=epn_cl_src.lc \
# 	 maketimecolumn=yes \
# 	 timecolumn=TIME \
# 	 timebinsize=100 \
# 	 makeratecolumn=yes
# dsplot table=epn_cl_src.lc x=TIME y=RATE

# Check for pile-up
epatplot set=epn_cl_src.evt \
	 plotfile=epn_pu.ps \
	 useplotfile=yes
# epatplot set=epn_cl_bkg.evt \
# 	 plotfile=epn_pu_bkg.ps \
# 	 useplotfile=yes

# Barycentric correction
cp epn_cl_src.evt epn_cl_src_bc.evt
barycen table=epn_cl_src_bc.evt:EVENTS \
        ephemeris=DE405 \
	srcra=${ra} \
	srcdec=${de}
# cp epn_cl_bkg.evt epn_cl_bkg_bc.evt
# barycen table=epn_cl_bkg_bc.evt:EVENTS \
#         ephemeris=DE405 \
# 	srcra=${ra} \
# 	srcdec=${de}

# PSF
# psfgen image=epn_cl.img \
#        energy=500 \
#        level=ELLBETA \
#        coordtype=POS \
#        x=${ra} \
#        y=${de} \
#        xsize=400 \
#        ysize=400 \
#        output=epn_0500.psf

# psfgen image=epn_cl.img \
#        energy=1000 \
#        level=ELLBETA \
#        coordtype=POS \
#        x=${ra} \
#        y=${de} \
#        xsize=400 \
#        ysize=400 \
#        output=epn_1000.psf

# psfgen image=epn_cl.img \
#        energy=2000 \
#        level=ELLBETA \
#        coordtype=POS \
#        x=${ra} \
#        y=${de} \
#        xsize=400 \
#        ysize=400 \
#        output=epn_2000.psf
