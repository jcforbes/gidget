#!/bin/bash
#SBATCH --job-name rf231
#SBATCH -n 10
#SBATCH -N 1-1
#SBATCH -t 2-00:00
#SBATCH -p shared 
#SBATCH --mem-per-cpu=1000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johncforbes@gmail.com

module load python/2.7.13-fasrc01
module load gsl

python exper.py --nproc 10 rf231
#mpirun -n 128 python mcmc_broad.py 
