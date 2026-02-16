#!/bin/bash
LOG_FILE=./sphere_global_losses_1x1_runTimeCompLocal.log
exec > >(tee ${LOG_FILE}) 2>&1
export JULIA_NUM_THREADS=1
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeCompLocal.jl --freq 100 --compute_full_solution true --mesh_file "sphere_1m_coarser"
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeCompLocal.jl --freq 1000 --compute_full_solution true --mesh_file "sphere_1m_coarser"
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeCompLocal.jl --freq 100 --compute_full_solution true --mesh_file "sphere_1m_coarse"
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeCompLocal.jl --freq 1000 --compute_full_solution true --mesh_file "sphere_1m_coarse"
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeCompLocal.jl --freq 100 --compute_full_solution true --mesh_file "sphere_1m"
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeCompLocal.jl --freq 1000 --compute_full_solution true --mesh_file "sphere_1m"
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeCompLocal.jl --freq 100 --compute_full_solution true --mesh_file "sphere_1m_fine"
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeCompLocal.jl --freq 1000 --compute_full_solution true --mesh_file "sphere_1m_fine"
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeCompLocal.jl --freq 100 --compute_full_solution true --mesh_file "sphere_1m_finer"
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeCompLocal.jl --freq 1000 --compute_full_solution true --mesh_file "sphere_1m_finer"
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeCompLocal.jl --freq 100 --compute_full_solution true --mesh_file "sphere_1m_4p5k"
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeCompLocal.jl --freq 1000 --compute_full_solution true --mesh_file "sphere_1m_4p5k"
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeCompLocal.jl --freq 100 --compute_full_solution true --mesh_file "sphere_1m_extremely_fine"
julia --project=.  ./examples/sphere_global_losses_1x1_runTimeCompLocal.jl --freq 1000 --compute_full_solution true --mesh_file "sphere_1m_extremely_fine"