#!/bin/bash -x

wrk_dir="/Users/silver/dat/nus/cow/"
ra="244.00093408"
de="22.26802508"
ids=(90401327002 90401327004 90401327006 90401327008)

################################################################

nid=${#ids[@]}
cd $wrk_dir

for id in $(seq 0 $(($nid-1)))
do
    cd ${ids[id]}
    echo ${ids[id]}
    
    # nupipeline indir=. steminputs=nu${ids[id]} outdir=./pro obsmode=SCIENCE saacalc=3 saamode=optimized tentacle=yes

    barycorr infile=./pro/nu${ids[id]}A01_cl.evt outfile=./pro/nu${ids[id]}A01_bc.evt orbitfiles=./auxil/nu${ids[id]}_orb.fits ra=${ra} dec=${de} refframe=ICRS
    barycorr infile=./pro/nu${ids[id]}B01_cl.evt outfile=./pro/nu${ids[id]}B01_bc.evt orbitfiles=./auxil/nu${ids[id]}_orb.fits ra=${ra} dec=${de} refframe=ICRS
    
# # THIS IS FOR SUBTRACTING BACKGROUND
#     nuproducts indir=./pro infile=./pro/nu${ids[id]}A01_bc.evt instrument=FPMA steminputs=nu${ids[id]} outdir=./pro srcregionfile=srcA.reg bkgregionfile=bkgA.reg rungrppha=yes grpmincounts=25 grppibadlow=35 grppibadhigh=1909 binsize=5808
#     nuproducts indir=./pro infile=./pro/nu${ids[id]}B01_bc.evt instrument=FPMB steminputs=nu${ids[id]} outdir=./pro srcregionfile=srcB.reg bkgregionfile=bkgB.reg rungrppha=yes grpmincounts=25 grppibadlow=35 grppibadhigh=1909 binsize=5808

#     /Applications/SAOImageDS9.app/Contents/MacOS/ds9 ./pro/nu${ids[id]}A01_bc.evt -scale linear -cmap Heat -regions ./srcA.reg -regions ./bkgA.reg -print destination file -print filename ./pro/nu${ids[id]}A01_img.ps -print -exit
#     /Applications/SAOImageDS9.app/Contents/MacOS/ds9 ./pro/nu${ids[id]}B01_bc.evt -scale linear -cmap Heat -regions ./srcB.reg -regions ./bkgB.reg -print destination file -print filename ./pro/nu${ids[id]}B01_img.ps -print -exit
    
#     for ff in ./pro/*.ps
#     do
#         ps2pdf ${ff}
#     done
#     rm ./pro/*.ps

    cd $wrk_dir
done

echo "EOF"
