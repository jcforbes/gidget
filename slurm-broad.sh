#!/bin/bash
#SBATCH --job-name brd130c
#SBATCH -n 128
#SBATCH -t 2-00:00
#SBATCH -p shared 
#SBATCH --mem-per-cpu=1000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johncforbes@gmail.com


module load gcc
module load openmpi/2.1.0-fasrc02
module load gsl/2.4-fasrc01
module load python/2.7.14-fasrc01

source activate ody


#python exper.py --nproc 16 rf264
mpirun -n 128 python mcmc_broad.py 
