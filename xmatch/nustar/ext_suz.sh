#!/bin/bash
# https://heasarc.gsfc.nasa.gov/docs/suzaku/aehp_data_analysis.html
# https://heasarc.gsfc.nasa.gov/docs/suzaku/analysis/abc/

# The filenames (except for some log files) use the following general convention:

# aeNNNNNNNNNiii_n_mmmmmmmm_ll.ext.gz
# where
# ae
# is short for Astro-E2, the initial name of Suzaku.
# NNNNNNNNN
# is the observation sequence number and is identical to the directory name.
# iii
# is the instrument specification. This string is set as follows: hxd=HXD, xi[0-3]=XIS-[0-3]. xis is used for files common to all the XIS units. This string can be omitted in files under the auxil and log directories.
# n
# ranges from 0 to 9 and indicates the RPT file number. The original telemetry file is divided into RPT files and more than one RPT can contribute to one observation. The value of 0 is used when the science file combines data from different RPT or if there is only one RPT file that contributes to that sequence. This number can be omitted in files under the auxil and log directories.
# mmmmmmmm
# is the file identifier. The string distinguishes between files from the same instrument.
# ll
# indicates the file level. For event files, the string can be ``uf'' or ``cl'' to indicate ``unfiltered'' or ``cleaned'' event files. It also can be ``bg,'' ``sk,'' ``sr,'' ``gso,'' ``pin,'' ``wel'' (products directory for both the XIS and HXD) or ``wam'' ( hk directory for the HXD). The string can be omitted.
# ext
# is the file extension. Currently it can take the values: ``evt'' (event files), ``gti'' (good time interval), ``hk'' (house keeping), ``ghf'' (gain history file), ``ght'' (gain history table), ``lc'' (light curve), ``pi'' (pulse invariant), ``html,'' ``log,'' ``com,'' ``att'' (attitude file), ``cat,'' ``ehk,'' ``mkf,'' ``orb,'' ``tim,'' ``img,'' and ``gif.''
# For more information on file names of the products of the pipeline processing, please refer to the documentation that can be found at
# http://heasarc.gsfc.nasa.gov/docs/suzaku/aehp_data_analysis.html.

# Aebarycen is a Suzaku specific FTOOL to perform barycentric corrections for Suzaku event and housekeeping files. See Terada et al., 2008, PASJ 60, S25 for a brief description of the barycentric correction tool (and a detailed description of the timing calibration of the HXD).
# Warning: In order for the barycentric correction tool to work correctly it is required to perform any good time interval filtering in xselect, i.e., removing times of telemetry saturation, before running aebarycen.


id="409049010"
wrk="/Users/silver/dat/suz/cx1"
out="/Users/silver/dat/suz/cx1/${id}_repro"
dat="/Users/silver/dat/suz/cx1/${id}"
ra="299.59031592"
de="35.20160611"

if [[ $out =~ " " ]]
then
    echo "Path to cwd should not contain spaces! (SAS doesn't allow this)"
    exit 1
fi
mkdir -p ${out}
cd ${wrk}

# aepipeline indir=${dat} \
#            outdir=${out} \
#            steminputs=ae${id} \
#            entry_stage=1 \
#            exit_stage=2 \
#            clobber=yes \
#            instrument=PIN

mv "${out}/ae409049010hxd_0_pinno_cl.evt" "${out}/ae409049010hxd_bc.evt"

aebarycen filelist="${out}/ae409049010hxd_bc.evt" \
          orbit="${out}/ae409049010.orb" \
          ra=${ra} \
          dec=${de}
