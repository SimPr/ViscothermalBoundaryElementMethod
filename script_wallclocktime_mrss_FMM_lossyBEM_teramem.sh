#!/bin/bash
#SBATCH -o ./myjob.%j.%N.out
#SBATCH -D ./
#SBATCH -J lossyBEM_wallclocktime_mrss_teramem 
#SBATCH --get-user-env
#SBATCH --clusters=inter
#SBATCH --partition=teramem_inter
#SBATCH --mem=1000000mb
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=simone.preuss@tum.de
#SBATCH --export=NONE
#SBATCH --time=24:00:00 
#--------------------------------------
module load slurm_setup
# Load other intended modules here....
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
/usr/bin/time -v julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp2c.jl --freq 100 --mesh_file "sphere_1m_coarser"
/usr/bin/time -v julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp2c.jl --freq 1000 --mesh_file "sphere_1m_coarser"
/usr/bin/time -v julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp2c.jl --freq 100 --mesh_file "sphere_1m_coarse"
/usr/bin/time -v julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp2c.jl --freq 1000 --mesh_file "sphere_1m_coarse"
/usr/bin/time -v julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp2c.jl --freq 100 --mesh_file "sphere_1m"
/usr/bin/time -v julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp2c.jl --freq 1000 --mesh_file "sphere_1m"
/usr/bin/time -v julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp2c.jl --freq 100 --mesh_file "sphere_1m_fine"
/usr/bin/time -v julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp2c.jl --freq 1000 --mesh_file "sphere_1m_fine"
/usr/bin/time -v julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp2c.jl --freq 100 --mesh_file "sphere_1m_finer"
/usr/bin/time -v julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp2c.jl --freq 1000 --mesh_file "sphere_1m_finer"
/usr/bin/time -v julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp2c.jl --freq 100 --mesh_file "sphere_1m_4p5k"
/usr/bin/time -v julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp2c.jl --freq 1000 --mesh_file "sphere_1m_4p5k"
/usr/bin/time -v julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp2c.jl --freq 100 --mesh_file "sphere_1m_extremely_fine"
/usr/bin/time -v julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp2c.jl --freq 1000 --mesh_file "sphere_1m_extremely_fine"
/usr/bin/time -v julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp2c.jl --freq 100 --mesh_file "sphere_1m_35k"
/usr/bin/time -v julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp2c.jl --freq 1000 --mesh_file "sphere_1m_35k"
/usr/bin/time -v julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp2c.jl --freq 100 --mesh_file "sphere_1m_77k"
/usr/bin/time -v julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp2c.jl --freq 1000 --mesh_file "sphere_1m_77k"