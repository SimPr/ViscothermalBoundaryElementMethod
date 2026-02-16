# Track time and memory for LGM setup (tLGM, equal for all), assembly of system matrix (tLGM_dense1x1 for 1n-direct, tLGM_dense4x4 for 4n-direct, tLGM_dense10x10 for 10n-direct)
# and solution (tIter for 1n-iter based on 1 GMRES iteration, tLGM_dense1x1 for 1n-direct, tLGM_dense4x4 for 4n-direct, tLGM_dense10x10 for 10n-direct).
# Track memory for storing LGM (varLGM, equal for all) and system matrix (varLGMdense1x1 for 1n-direct, varLGMdense4x4 for 4n-direct, varLGMdense10x10 for 10n-direct).
# Call with julia --project=.  ./examples/sphere_global_losses_MemComp.jl --freq 500 --mesh_file "sphere_1m_coarser"

using ArgParse
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--freq"
            help = "frequency"
            default = 1000.0
            arg_type = Float64
        "--mesh_file"
            help = "mesh file"
            default = "sphere_1m"
            arg_type = String

    end
    return parse_args(s)
end
function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
    return parsed_args
end

intputArguments = main()
freq = intputArguments["freq"]
mesh_file = intputArguments["mesh_file"]

#==========================================================================================
                            Adding Related Packages
==========================================================================================#
using LinearAlgebra
using BoundaryIntegralEquations
using IterativeSolvers
using BenchmarkTools
using JLD2
using Statistics, DataFrames


### Auxiliary functions for benchmarking
preprocess_trial(t::BenchmarkTools.Trial, id::AbstractString) = (id=id,
        minimum=minimum(t.times),
        median=median(t.times),
        maximum=maximum(t.times),
        allocations=t.allocs,
        memory_estimate=t.memory)

#=============================ß============================================================
                                Loading Mesh
==========================================================================================#
geometry_orders     = [:linear,:quadratic];
tri_physics_orders  = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic];
# Triangular Meshes

mesh_path = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes");

#for mesh_file in mesh_files
tri_mesh_file = joinpath(mesh_path,mesh_file);
mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[2],
                                                    physics_order=tri_physics_orders[2]);
#==========================================================================================
                            Creating excitation vector
==========================================================================================#
xyzb = mesh.sources;
M  = size(xyzb,2);
u₀ = 1e-2;
v0 = [zeros(2M); u₀*ones(M)];
#===========================================================================================
                        BEM matrix assembly and (iterative) solution of the 1-variable system
===========================================================================================#

output = DataFrame()
@info "Computing LGM"
push!(output, preprocess_trial(@benchmark(LossyGlobalOuter(mesh,freq;fmm_on=false,depth=1,n=3,progress=false),evals=5), "tLGM"))
LGM = LossyGlobalOuter(mesh,freq;fmm_on=false,depth=1,n=3,progress=false);

info = Base.summarysize(LGM);

# 1x1 iterative format
rhs = LGM.Ga*gmres(LGM.inner,(LGM.Dr*v0 - LGM.Nd'*gmres(LGM.Gv,LGM.Hv*v0;verbose=false));verbose=true);
@info "Assembling and solving iterative system"
push!(output, preprocess_trial(@benchmark(gmres(LGM,rhs;verbose=false,maxiter=1),evals=5), "tIter"))

# 1x1 dense
@info "Assembling dense 1x1 system"
push!(output, preprocess_trial(@benchmark(BoundaryIntegralEquations._full1_new(LGM),evals=5), "tLGM_dense1x1"))
LGM_dense  = BoundaryIntegralEquations._full1_new(LGM);

info1 = Base.summarysize(LGM_dense);

@info "Solving dense 1x1 system"
push!(output, preprocess_trial(@benchmark(LGM_dense\rhs,evals=5), "tDense1x1"))


# 4x4 dense
rhs4 = [zeros(M); v0];
@info "Assembling dense 4x4 system"
push!(output, preprocess_trial(@benchmark(BoundaryIntegralEquations._full4(LGM),evals=5), "tLGM_dense4x4"))
LGM_dense4  = BoundaryIntegralEquations._full4(LGM);

info4 = Base.summarysize(LGM_dense4);

@info "Solving dense 4x4 system"
push!(output, preprocess_trial(@benchmark(LGM_dense4\rhs4,evals=5), "tDense4x4"))

# 10x10 dense
rhs10 = [zeros(7M);v0];
@info "Assembling dense 10x10 system"
push!(output, preprocess_trial(@benchmark(BoundaryIntegralEquations._full10(LGM),evals=5), "tLGM_dense10x10"))
LGM_dense10  = BoundaryIntegralEquations._full10(LGM);

info10 = Base.summarysize(LGM_dense10);

@info "Solving dense 10x10 system"
push!(output, preprocess_trial(@benchmark(LGM_dense10\rhs10,evals=5), "tDense10x10"))


jldsave("memory_$(M)DOFs_$(Int(freq))Hz.JLD2", 
memory=output,
varLGM=info,
varLGMdense1x1=info1,
varLGMdense4x4=info4,
varLGMdense10x10=info10,
)

println(output)

println("varInfo tLGM, +tLGM_dense1x1, +tLGM_dense4x4, +tLGM_dense10x10")
println(info)
println(info1)
println(info4)
println(info10)