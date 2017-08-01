#!/bin/bash
#SBATCH --job-name brd21c
#SBATCH -n 256
#SBATCH -t 4-00:00
#SBATCH -p general
#SBATCH --mem-per-cpu=1000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johncforbes@gmail.com

module load python/2.7.13-fasrc01

mpirun -n 256 python mcmc_broad.py 
