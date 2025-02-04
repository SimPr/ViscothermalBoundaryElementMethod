#!/bin/bash
#SBATCH -o ./myjob.%j.%N.out
#SBATCH -D ./
#SBATCH -J lossyIGABEM_results_teramem
#SBATCH --get-user-env
#SBATCH --clusters=inter
#SBATCH --partition=teramem_inter
#SBATCH --mem=300000mb
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=simone.preuss@tum.de
#SBATCH --export=NONE
#SBATCH --time=24:00:00 
#--------------------------------------
module load slurm_setup
# Load other intended modules here....
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 100 --compute_full_solution true --mesh_file "sphere_1m_coarser"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 100 --compute_full_solution true --mesh_file "sphere_1m_coarse"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 100 --compute_full_solution true --mesh_file "sphere_1m"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 100 --compute_full_solution true --mesh_file "sphere_1m_fine"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 100 --compute_full_solution true --mesh_file "sphere_1m_fine2"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 200 --compute_full_solution true --mesh_file "sphere_1m_coarser"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 200 --compute_full_solution true --mesh_file "sphere_1m_coarse"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 200 --compute_full_solution true --mesh_file "sphere_1m"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 200 --compute_full_solution true --mesh_file "sphere_1m_fine"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 200 --compute_full_solution true --mesh_file "sphere_1m_fine2"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 300 --compute_full_solution true --mesh_file "sphere_1m_coarser"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 300 --compute_full_solution true --mesh_file "sphere_1m_coarse"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 300 --compute_full_solution true --mesh_file "sphere_1m"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 300 --compute_full_solution true --mesh_file "sphere_1m_fine"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 300 --compute_full_solution true --mesh_file "sphere_1m_fine2"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 400 --compute_full_solution true --mesh_file "sphere_1m_coarser"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 400 --compute_full_solution true --mesh_file "sphere_1m_coarse"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 400 --compute_full_solution true --mesh_file "sphere_1m"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 400 --compute_full_solution true --mesh_file "sphere_1m_fine"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 400 --compute_full_solution true --mesh_file "sphere_1m_fine2"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 500 --compute_full_solution true --mesh_file "sphere_1m_coarser"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 500 --compute_full_solution true --mesh_file "sphere_1m_coarse"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 500 --compute_full_solution true --mesh_file "sphere_1m"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 500 --compute_full_solution true --mesh_file "sphere_1m_fine"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 500 --compute_full_solution true --mesh_file "sphere_1m_fine2"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 600 --compute_full_solution true --mesh_file "sphere_1m_coarser"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 600 --compute_full_solution true --mesh_file "sphere_1m_coarse"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 600 --compute_full_solution true --mesh_file "sphere_1m"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 600 --compute_full_solution true --mesh_file "sphere_1m_fine"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 600 --compute_full_solution true --mesh_file "sphere_1m_fine2"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 700 --compute_full_solution true --mesh_file "sphere_1m_coarser"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 700 --compute_full_solution true --mesh_file "sphere_1m_coarse"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 700 --compute_full_solution true --mesh_file "sphere_1m"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 700 --compute_full_solution true --mesh_file "sphere_1m_fine"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 700 --compute_full_solution true --mesh_file "sphere_1m_fine2"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 800 --compute_full_solution true --mesh_file "sphere_1m_coarser"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 800 --compute_full_solution true --mesh_file "sphere_1m_coarse"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 800 --compute_full_solution true --mesh_file "sphere_1m"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 800 --compute_full_solution true --mesh_file "sphere_1m_fine"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 800 --compute_full_solution true --mesh_file "sphere_1m_fine2"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 900 --compute_full_solution true --mesh_file "sphere_1m_coarser"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 900 --compute_full_solution true --mesh_file "sphere_1m_coarse"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 900 --compute_full_solution true --mesh_file "sphere_1m"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 900 --compute_full_solution true --mesh_file "sphere_1m_fine"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 900 --compute_full_solution true --mesh_file "sphere_1m_fine2"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 1000 --compute_full_solution true --mesh_file "sphere_1m_coarser"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 1000 --compute_full_solution true --mesh_file "sphere_1m_coarse"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 1000 --compute_full_solution true --mesh_file "sphere_1m"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 1000 --compute_full_solution true --mesh_file "sphere_1m_fine"
julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 1000 --compute_full_solution true --mesh_file "sphere_1m_fine2"