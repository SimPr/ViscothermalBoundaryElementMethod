#!/bin/bash
#SBATCH -o ./myjob.%j.%N.out
#SBATCH -D ./
#SBATCH -J lossyFMBEM_resultsFMM_teramem
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
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 100 --mesh_file "sphere_1m_coarser" --porder 2
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 1000 --mesh_file "sphere_1m_coarser" --porder 2
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 100 --mesh_file "sphere_1m_coarse" --porder 2
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 1000 --mesh_file "sphere_1m_coarse" --porder 2
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 100 --mesh_file "sphere_1m" --porder 2
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 1000 --mesh_file "sphere_1m" --porder 2
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 100 --mesh_file "sphere_1m_fine" --porder 2
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 1000 --mesh_file "sphere_1m_fine" --porder 2
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 100 --mesh_file "sphere_1m_finer" --porder 2
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 1000 --mesh_file "sphere_1m_finer" --porder 2
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 100 --mesh_file "sphere_1m_4p5k" --porder 2
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 1000 --mesh_file "sphere_1m_4p5k" --porder 2
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 100 --mesh_file "sphere_1m_extremely_fine" --porder 2
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 1000 --mesh_file "sphere_1m_extremely_fine" --porder 2
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 100 --mesh_file "sphere_1m_35k" --porder 2
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 1000 --mesh_file "sphere_1m_35k" --porder 2
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 100 --mesh_file "sphere_1m_77k" --porder 2
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 1000 --mesh_file "sphere_1m_77k" --porder 2
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 100 --mesh_file "sphere_1m_coarser" --porder 1
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 1000 --mesh_file "sphere_1m_coarser" --porder 1
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 100 --mesh_file "sphere_1m_coarse" --porder 1
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 1000 --mesh_file "sphere_1m_coarse" --porder 1
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 100 --mesh_file "sphere_1m" --porder 1
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 1000 --mesh_file "sphere_1m" --porder 1
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 100 --mesh_file "sphere_1m_fine" --porder 1
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 1000 --mesh_file "sphere_1m_fine" --porder 1
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 100 --mesh_file "sphere_1m_finer" --porder 1
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 1000 --mesh_file "sphere_1m_finer" --porder 1
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 100 --mesh_file "sphere_1m_4p5k" --porder 1
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 1000 --mesh_file "sphere_1m_4p5k" --porder 1
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 100 --mesh_file "sphere_1m_extremely_fine" --porder 1
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 1000 --mesh_file "sphere_1m_extremely_fine" --porder 1
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 100 --mesh_file "sphere_1m_35k" --porder 1
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 1000 --mesh_file "sphere_1m_35k" --porder 1
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 100 --mesh_file "sphere_1m_77k" --porder 1
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 1000 --mesh_file "sphere_1m_77k" --porder 1
