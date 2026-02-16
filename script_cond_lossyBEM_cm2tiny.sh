#!/bin/bash
#SBATCH -o ./myjob.%j.%N.out
#SBATCH -D ./
# --- By using relative path names the script needs to be submitted from the work directory ---
#SBATCH -J lossyBEM_benchmark_cm2tiny
#SBATCH --cluster=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
# --- multiples of 28 for CM2 cluster with Intel Haswell processors ---
#SBATCH --mail-type=ALL
#SBATCH --mail-user=simone.preuss@tum.de
#SBATCH --time=12:00:00
#----------------------------------------------------
module load slurm_setup
# Load other intended modules here....
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
julia --project=.  ./examples/paper1_comparisons1.jl