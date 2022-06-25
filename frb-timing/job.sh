#!/bin/bash -x

# job name
#SBATCH -J frb
# account
#SBATCH -A 2018-101
# email notification
#SBATCH --mail-user=dalp@kth.se
#SBATCH --mail-type=ALL
# wall-clock time to be given to this job
#SBATCH -t 05:00:00
# Number of nodes
#SBATCH --nodes=1
# set tasks per node to 24 in order to disable hyperthreading
#SBATCH --ntasks-per-node=24
#SBATCH -e err.log
#SBATCH -o out.log
#SBATCH --mem=1000000

module add i-compilers intelmpi
./compile.sh clean
./compile.sh tegner

ff=($(ls ./dat/zhang18b.txt))
nf=${#ff[*]}
rm -rf out*

for id in $(seq 0 $(($nf-1)))
do
    INP_FIL=${ff[${id}]}
    mkdir -p out
    mpirun -np 24 ./bin/tim $INP_FIL $(wc -l < $INP_FIL) > out.log 2>&1
    mv out/dis.txt out/dis_${id}.txt
    mv out/top.txt out/top_${id}.txt
done
