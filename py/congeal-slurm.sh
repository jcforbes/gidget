#!/bin/bash
#SBATCH --job-name congeal
#SBATCH -n 1
#SBATCH -t 2-00:00
#SBATCH -p shared
#SBATCH --mem-per-cpu=10000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johncforbes@gmail.com

module unload python/2.7.6-fasrc01

module load python/2.7.13-fasrc01

cd $GIDGETDIR/py
python congeal.py > estd.${SLURM_JOB_ID}.out 

