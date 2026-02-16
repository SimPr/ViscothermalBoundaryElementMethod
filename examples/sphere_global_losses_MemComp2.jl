# Track time and memory for solution of 1n-iter based on 20 GMRES iterations as restart happens at 20 (tIter) (only relevant for 1000Hz as 100Hz study requires less than 20 iterations)
# Call with julia --project=.  ./examples/sphere_global_losses_MemComp2.jl --freq 500 --mesh_file "sphere_1m_coarser"

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
LGM = LossyGlobalOuter(mesh,freq;fmm_on=false,depth=1,n=3,progress=false);


# 1x1 iterative format
rhs = LGM.Ga*gmres(LGM.inner,(LGM.Dr*v0 - LGM.Nd'*gmres(LGM.Gv,LGM.Hv*v0;verbose=false));verbose=false);
@info "Assembling and solving iterative system"
push!(output, preprocess_trial(@benchmark(gmres(LGM,rhs;verbose=false,maxiter=20),evals=5), "tIter"))


jldsave("memoryMaxIter_$(M)DOFs_$(Int(freq))Hz.JLD2", 
memory=output,
)

println(output)