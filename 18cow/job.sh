#!/bin/bash -l

# job name
#SBATCH -J cow
# account
#SBATCH -A 2018-101
# email notification
#SBATCH --mail-user=dalp@kth.se
#SBATCH --mail-type=ALL
# wall-clock time to be given to this job
#SBATCH -t 75:00:00
# Number of nodes
#SBATCH --nodes=2
# set tasks per node to 24 in order to disable hyperthreading
#SBATCH --ntasks-per-node=24
#SBATCH -e err.log
#SBATCH -o out.log

module add i-compilers intelmpi
./compile.sh clean
./compile.sh tegner

ff=($(ls ./dat/12+_*))
nf=21
rm out/*

for id in $(seq 0 $(($nf-1)))
do
    INP_FIL=${ff[${id}]}
    mkdir -p out
    mpirun -np 48 ./bin/cow $INP_FIL $(wc -l < $INP_FIL)
    id2=$(printf %03d ${id})
    mv out/dis.txt out/dis_${id2}.txt
    mv out/PP.bin out/PP_${id2}.bin
    mv out/Pd.bin out/Pd_${id2}.bin
    mv out/HH.bin out/HH_${id2}.bin
done
