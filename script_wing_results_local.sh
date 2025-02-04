#!/bin/bash
LOG_FILE=./wing_global_losses_1x1_runTimeCompLocal.log
exec > >(tee ${LOG_FILE}) 2>&1
export JULIA_NUM_THREADS=6
/usr/bin/time -v julia --project=.  ./examples/wing_global_losses_1x1_exteriorProblem.jl --freq 500 --compute_full_solution true --mesh_file "wing_extremely_fine"