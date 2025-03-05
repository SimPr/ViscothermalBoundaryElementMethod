#!/bin/bash
#SBATCH -o ./myjob.%j.%N.out
#SBATCH -D ./
#SBATCH -J lossyIGABEM_FreqStudy_cm2_tiny
#SBATCH --cluster=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
# --- multiples of 28 for CM2 cluster with Intel Haswell processors ---
#SBATCH --mail-type=ALL
#SBATCH --mail-user=simone.preuss@tum.de
#SBATCH --time=10:00:00
#--------------------------------------
module load slurm_setup
# Load other intended modules here....
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
END=1000
for ((i=100;i<=END;i++)); do
    julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq $i --compute_full_solution true --mesh_file "sphere_1m_fine"
done
