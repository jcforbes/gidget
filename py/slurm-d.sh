#!/bin/bash
#SBATCH --job-name dyn02a
#SBATCH -n 15
#SBATCH -N 1
#SBATCH -t 4-12:00
#SBATCH -p shared
#SBATCH --mem=20000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johncforbes@gmail.com



module load gcc
module load openmpi/2.1.0-fasrc02
module load gsl/2.4-fasrc01
module load python/2.7.14-fasrc01

source activate ody

cd $GIDGETDIR/py
#mpirun -n 15 python broad_svm.py --runEmcee --seedWith='intq02951_fakemcmc140a_restart.pickle' --fakemcmcName='fakemcmc140b' --logMh0=10.0 > estd.${SLURM_JOB_ID}.out 

python broad_svm.py --runDynesty --fakemcmcName='dyn02' --logMh0=10.0 > estd.${SLURM_JOB_ID}.out 
