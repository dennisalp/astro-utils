#!/bin/bash

wrk_dir=/Volumes/pow/dat/nus/14j/

################################################################

cd $wrk_dir
ids=($(ls -d */ | sed 's#/##'))
nid=${#ids[*]}

for id in $(seq 2 $(($nid-4)))
#for id in $(seq 0 $(($nid-1)))
#for id in $(seq 0 $(($nid-7)))
do
    cd ${ids[id]}
    echo ${ids[id]}
    
#    nupipeline indir=. steminputs=nu${ids[id]} outdir=./pro obsmode=SCIENCE

# THIS IS FOR SUBTRACTING BACKGROUND
#    nuproducts indir=./pro infile=./pro/nu${ids[id]}A01_cl.evt instrument=FPMA steminputs=nu${ids[id]} outdir=./pro srcregionfile=srcA.reg bkgregionfile=bkgA.reg rungrppha=yes grpmincounts=25 grppibadlow=35 grppibadhigh=1909 binsize=5808 pilow=1460
#    nuproducts indir=./pro infile=./pro/nu${ids[id]}B01_cl.evt instrument=FPMB steminputs=nu${ids[id]} outdir=./pro srcregionfile=srcB.reg bkgregionfile=bkgB.reg rungrppha=yes grpmincounts=25 grppibadlow=35 grppibadhigh=1909 binsize=5808 pilow=1460

    # THIS IS FOR FITTING BACKGROUND
    # First, extract the source
    nuproducts indir=./pro infile=./pro/nu${ids[id]}A01_cl.evt instrument=FPMA steminputs=nu${ids[id]} outdir=./pro srcregionfile=nu${ids[id]}A01_src.reg rungrppha=yes grpmincounts=25 grppibadlow=35 grppibadhigh=1909 binsize=512 pilow=1460 imagefile=NONE lcfile=NONE bkgextract=no
    nuproducts indir=./pro infile=./pro/nu${ids[id]}B01_cl.evt instrument=FPMB steminputs=nu${ids[id]} outdir=./pro srcregionfile=nu${ids[id]}B01_src.reg rungrppha=yes grpmincounts=25 grppibadlow=35 grppibadhigh=1909 binsize=512 pilow=1460 imagefile=NONE lcfile=NONE bkgextract=no

    # Second, extract the background
    nuproducts indir=./pro infile=./pro/nu${ids[id]}A01_cl.evt instrument=FPMA steminputs=nu${ids[id]} outdir=./pro srcregionfile=nu${ids[id]}A01_bkg.reg rungrppha=yes grpmincounts=25 grppibadlow=35 grppibadhigh=1909 binsize=512 pilow=1460 imagefile=NONE lcfile=NONE bkgextract=no extended=yes boxsize=20 phafile=nu${ids[id]}A01_bk.pha outarffile=nu${ids[id]}A01_bk.arf outrmffile=nu${ids[id]}A01_bk.rmf grpphafile=nu${ids[id]}A01_bk_grp.pha vignflag=no apstopflag=no detabsflag=no psfflag=no grflag=no
    nuproducts indir=./pro infile=./pro/nu${ids[id]}B01_cl.evt instrument=FPMB steminputs=nu${ids[id]} outdir=./pro srcregionfile=nu${ids[id]}B01_bkg.reg rungrppha=yes grpmincounts=25 grppibadlow=35 grppibadhigh=1909 binsize=512 pilow=1460 imagefile=NONE lcfile=NONE bkgextract=no extended=yes boxsize=20 phafile=nu${ids[id]}B01_bk.pha outarffile=nu${ids[id]}B01_bk.arf outrmffile=nu${ids[id]}B01_bk.rmf grpphafile=nu${ids[id]}B01_bk_grp.pha vignflag=no apstopflag=no detabsflag=no psfflag=no grflag=no
    
    /Applications/SAOImage\ DS9.app/Contents/MacOS/ds9 ./pro/nu${ids[id]}A01_cl.evt -scale log -cmap Heat -regions ./nu${ids[id]}A01_src.reg -regions ./nu${ids[id]}A01_bkg.reg -print destination file -print filename ./pro/nu${ids[id]}A01_img.ps -print -exit
    /Applications/SAOImage\ DS9.app/Contents/MacOS/ds9 ./pro/nu${ids[id]}B01_cl.evt -scale log -cmap Heat -regions ./nu${ids[id]}B01_src.reg -regions ./nu${ids[id]}B01_bkg.reg -print destination file -print filename ./pro/nu${ids[id]}B01_img.ps -print -exit
    
    for ff in ./pro/*.ps
    do
        ps2pdf ${ff}
    done
    rm ./pro/*.ps

    cd $wrk_dir
done

echo "EOF"
