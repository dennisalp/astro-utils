# ssh stuff
scp ./out_light/* dalp@t04n28.pdc.kth.se:/cfs/klemming/scratch/d/dalp/out_light/
cd /cfs/klemming/scratch/d/dalp/
ssh dalp@tegner.pdc.kth.se
kinit -f -l 5h30m dalp@NADA.KTH.SE

# running python on interactive node on tegner
ssh t02n13.pdc.kth.se
module load anaconda/py35/4.2.0
source activate custom

module load anaconda/py27/5.0.1
python some_script.py

# slurm on pdc
salloc -A 2018-101 -t 00:30:00 -N 1 --mem=2000000
mpirun -np 48 ./program
srun -n 1 ./program
aprun -n 1 python some_script.py
echo $SLURM_NODELIST
projinfo
sbatch job.sh
