#!/bin/bash -l

# job name
#SBATCH -J vis
# account
#SBATCH -A 2018-101
# email notification
#SBATCH --mail-user=dalp@kth.se
#SBATCH --mail-type=ALL
# wall-clock time to be given to this job
#SBATCH -t 15:00:00
# Number of nodes
#SBATCH --nodes=1
# set tasks per node to 24 in order to disable hyperthreading
#SBATCH --ntasks-per-node=1
#SBATCH -n 1
#SBATCH -e err.log
#SBATCH -o out.log

module load anaconda/py35/4.2.0
source activate custom

python src/vis_serial.py
