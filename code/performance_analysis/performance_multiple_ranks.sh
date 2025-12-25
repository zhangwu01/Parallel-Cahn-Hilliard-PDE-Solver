#!/usr/bin/env bash
# File       : performance_multiple_ranks.sh
# Job script to measure performance for MPI multi process performance
# Copyright 2023 Harvard University. All Rights Reserved.
#SBATCH --job-name=cs205_team13_check_ranks_anqi
#SBATCH --output=team13_check%j_ranks_anqi.out
#SBATCH --mem-per-cpu=2048
#SBATCH --nodes=1
#SBATCH --tasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=00:45:00
#SBATCH --mail-user=anqichen@g.harvard.edu
#SBATCH --mail-type=END

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



#threads=(1 1 1 4 4 4 8 8 8 12 12 12 20 20 20 28 28 28)
threads=(1)
#ranks=(1)
ranks=(1 1 1 4 4 4 9 9 9 16 16 16)


for rank_count in "${ranks[@]}"
do
    for thread_count in "${threads[@]}"
    do
        export OMP_NUM_THREADS=$thread_count
        echo "Program started with $rank_count ranks, $thread_count threads..."
        time srun -n$rank_count ./main 100
        time srun -n$rank_count ./main 350
        echo "Program completed with $rank_count ranks, $thread_count threads."
    done
   
done

# clean up
make clean
