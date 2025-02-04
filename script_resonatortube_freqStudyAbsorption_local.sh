#!/bin/bash
LOG_FILE=./resonatortube_global_losses_1x1_FreqStudyAbsorption.log
exec > >(tee ${LOG_FILE}) 2>&1
export JULIA_NUM_THREADS=6
END=301
for ((i=199;i<=END;i++)); do
    julia --project=.  ./examples/resonatortube_global_losses+lossless_1x1_absorp.jl --freq $i --mesh_file "resonator_tube"
done