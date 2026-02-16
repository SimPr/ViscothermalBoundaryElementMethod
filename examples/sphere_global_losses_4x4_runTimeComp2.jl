# Compute solution for pa (4n-direct) without extra steps to track total runtime and memory
# Call with /usr/bin/time -v julia --project=.  ./examples/sphere_global_losses_4x4_runTimeComp2.jl --freq 100 --mesh_file "sphere_1m_coarser"

using ArgParse
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--freq"
            help = "frequency"
            default = 100.0
            arg_type = Float64
        "--mesh_file"
            help = "mesh file"
            default = "sphere_1m_coarser"
            arg_type = String
    end
    return parse_args(s)
end
function main()
    return parse_commandline()
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
#=============================ß============================================================
                                Loading Mesh
==========================================================================================#
geometry_orders     = [:linear,:quadratic];
tri_physics_orders  = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic];
mesh_path = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes");
tri_mesh_file = joinpath(mesh_path,mesh_file);
mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[2],
                                                    physics_order=tri_physics_orders[2]);
#==========================================================================================
                            Creating excitation vector
==========================================================================================#
M  = size(mesh.sources,2);
u₀ = 1e-2;
v0 = [zeros(2M); u₀*ones(M)];
#===========================================================================================
                        BEM matrix assembly and (iterative) solution of the 1-variable system
===========================================================================================#
@info "Computing setup pa 4x4 dense $(M)DOFs freq$(Int(freq))"
LGM = LossyGlobalOuter(mesh,freq;fmm_on=false,depth=1,n=3,progress=false);
rhs = [zeros(M); v0];
LGM_dense  = BoundaryIntegralEquations._full4(LGM);
sol = LGM_dense\rhs;