#!/bin/bash
#SBATCH -o ./myjob.%j.%N.out
#SBATCH -D ./
#SBATCH -J lossyBEM_memory_teramem 
#SBATCH --get-user-env
#SBATCH --clusters=inter
#SBATCH --partition=teramem_inter
#SBATCH --mem=1000000mb
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=simone.preuss@tum.de
#SBATCH --export=NONE
#SBATCH --time=96:00:00 
#--------------------------------------
module load slurm_setup
# Load other intended modules here....
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
julia --project=.  ./examples/sphere_global_losses_MemComp.jl --freq 100 --mesh_file "sphere_1m_coarser"
julia --project=.  ./examples/sphere_global_losses_MemComp.jl --freq 100 --mesh_file "sphere_1m_coarse"
julia --project=.  ./examples/sphere_global_losses_MemComp.jl --freq 100 --mesh_file "sphere_1m"
julia --project=.  ./examples/sphere_global_losses_MemComp.jl --freq 100 --mesh_file "sphere_1m_finer"
julia --project=.  ./examples/sphere_global_losses_MemComp.jl --freq 100 --mesh_file "sphere_1m_fine"
julia --project=.  ./examples/sphere_global_losses_MemComp.jl --freq 100 --mesh_file "sphere_1m_4p5k"
julia --project=.  ./examples/sphere_global_losses_MemComp.jl --freq 100 --mesh_file "sphere_1m_extremely_fine"
julia --project=.  ./examples/sphere_global_losses_MemComp.jl --freq 1000 --mesh_file "sphere_1m_coarser"
julia --project=.  ./examples/sphere_global_losses_MemComp.jl --freq 1000 --mesh_file "sphere_1m_coarse"
julia --project=.  ./examples/sphere_global_losses_MemComp.jl --freq 1000 --mesh_file "sphere_1m"
julia --project=.  ./examples/sphere_global_losses_MemComp.jl --freq 1000 --mesh_file "sphere_1m_finer"
julia --project=.  ./examples/sphere_global_losses_MemComp.jl --freq 1000 --mesh_file "sphere_1m_fine"
julia --project=.  ./examples/sphere_global_losses_MemComp.jl --freq 1000 --mesh_file "sphere_1m_4p5k"
julia --project=.  ./examples/sphere_global_losses_MemComp.jl --freq 1000 --mesh_file "sphere_1m_extremely_fine"
julia --project=.  ./examples/sphere_global_losses_MemComp2.jl --freq 1000 --mesh_file "sphere_1m_coarser"
julia --project=.  ./examples/sphere_global_losses_MemComp2.jl --freq 1000 --mesh_file "sphere_1m_coarse"
julia --project=.  ./examples/sphere_global_losses_MemComp2.jl --freq 1000 --mesh_file "sphere_1m"
julia --project=.  ./examples/sphere_global_losses_MemComp2.jl --freq 1000 --mesh_file "sphere_1m_finer"
julia --project=.  ./examples/sphere_global_losses_MemComp2.jl --freq 1000 --mesh_file "sphere_1m_fine"
julia --project=.  ./examples/sphere_global_losses_MemComp2.jl --freq 1000 --mesh_file "sphere_1m_4p5k"
julia --project=.  ./examples/sphere_global_losses_MemComp2.jl --freq 1000 --mesh_file "sphere_1m_extremely_fine"