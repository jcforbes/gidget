#!/bin/bash
#SBATCH --job-name fmc110b
#SBATCH -n 15 
#SBATCH -t 1-06:00
#SBATCH -p shared
#SBATCH --mem-per-cpu=10000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johncforbes@gmail.com



module load gcc
module load openmpi/2.1.0-fasrc02
module load gsl/2.4-fasrc01
module load python/2.7.14-fasrc01

source activate ody

cd $GIDGETDIR/py
mpirun -n 15 python broad_svm.py --runEmcee --seedWith='' --fakemcmcName='fakemcmc110b' --logMh0=10.0 > estd.${SLURM_JOB_ID}.out 

