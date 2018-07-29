#!/bin/bash
#SBATCH --job-name fmc162a
#SBATCH -n 20
#SBATCH -t 2-06:00
#SBATCH -p shared
#SBATCH --mem-per-cpu=15000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johncforbes@gmail.com



module load gcc
module load openmpi/2.1.0-fasrc02
module load gsl/2.4-fasrc01
module load python/2.7.14-fasrc01

source activate ody

cd $GIDGETDIR/py
mpirun -n 20 python broad_svm.py --runEmcee --seedWith='intq03501_fakemcmc161b_restart.pickle' --fakemcmcName='fakemcmc162a' --logMh0Min=11.0 --logMh0Max=12.0 > estd.${SLURM_JOB_ID}.out 
#mpirun -n 15 python broad_svm.py --runEmcee --fakemcmcName='fakemcmc152a' --logMh0=12.0 > estd.${SLURM_JOB_ID}.out 

