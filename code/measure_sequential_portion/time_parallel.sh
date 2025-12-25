#!/usr/bin/env bash
# File       : submit_check.sh
# Description: Submit a test job to check your implementation
# Copyright 2023 Harvard University. All Rights Reserved.
#SBATCH --job-name=cs205_team13_check
#SBATCH --output=team13_check%j.out
#SBATCH --mem-per-cpu=2048
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --mail-user=anqichen@g.harvard.edu
#SBATCH --mail-type=END


module purge
module load python/3.9.12-fasrc01 gcc/12.1.0-fasrc01 openmpi/4.1.3-fasrc01

# build
root=$(pwd -P)
cwd=team13_check_${SLURM_JOBID}
mkdir -p ${cwd}
cd ${cwd}
cp -t . ${root}/*.h ${root}/*.cpp ${root}/*.py ${root}/Makefile
make clean_pedantic main

# run
export OMP_DYNAMIC='false'
export OMP_PLACE='core'
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
time srun ./main

# create figure
python plot.py
make clean
