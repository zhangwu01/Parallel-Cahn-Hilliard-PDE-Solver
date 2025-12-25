#!/usr/bin/env bash
# File       : example_job.sh
# Description: An example job script for a test case
# Copyright 2023 Harvard University. All Rights Reserved.
#SBATCH --job-name=cs205_team13_check
#SBATCH --output=team13_check%j.out
#SBATCH --mem-per-cpu=2048
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --time=00:15:00

module purge
module load python/3.9.12-fasrc01 gcc/12.1.0-fasrc01 openmpi/4.1.3-fasrc01

# build
root=$(pwd -P)
cwd=team13_check_${SLURM_JOBID}
mkdir -p ${cwd}
cd ${cwd}
mkdir -p steps
cp -t . ${root}/*.h ${root}/*.cpp ${root}/*.py ${root}/Makefile
make clean_pedantic main

# run
export OMP_DYNAMIC='false'
export OMP_PLACE='core'
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
time srun ./main 50 # command line argument for grid size

# create figure
python plot.py
make clean