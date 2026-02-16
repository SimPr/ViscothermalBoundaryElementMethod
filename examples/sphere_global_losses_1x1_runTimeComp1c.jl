# Save pa plot and solution for pa (1n-iter) and paFMM(1n-FMM) + analytical reference
# Call with julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1c.jl --freq 100 --mesh_file "sphere_1m_coarser" --porder 1

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
        "--porder"
            help = "order of physics discretization"
            default = 2
            arg_type = Int

    end
    return parse_args(s)
end
function main()
    return parse_commandline()
end

intputArguments = main()
freq = intputArguments["freq"]
mesh_file = intputArguments["mesh_file"]
porder = intputArguments["porder"]
#==========================================================================================
                            Adding Related Packages
==========================================================================================#
using LinearAlgebra
using BoundaryIntegralEquations
using IterativeSolvers
using JLD2
using Plots
#=============================ß============================================================
                                Loading Mesh
==========================================================================================#
geometry_orders     = [:linear,:quadratic];
tri_physics_orders  = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic];
mesh_path = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes");
tri_mesh_file = joinpath(mesh_path,mesh_file);
mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[2],
                                                    physics_order=tri_physics_orders[porder]);
#==========================================================================================
                            Creating excitation vector
==========================================================================================#
M  = size(mesh.sources,2);
u₀ = 1e-2;
v0 = [zeros(2M); u₀*ones(M)];
#===========================================================================================
                        FMM/full BEM matrix assembly and iterative solution of the 1-variable system
===========================================================================================#
@info "Computing setup pa 1x1 iterative $(M)DOFs freq$(Int(freq))"
LGM = LossyGlobalOuter(mesh,freq;fmm_on=false,depth=1,n=3,progress=false);
rhs = LGM.Ga*gmres(LGM.inner,(LGM.Dr*v0 - LGM.Nd'*gmres(LGM.Gv,LGM.Hv*v0;verbose=false));verbose=false);
pa = gmres(LGM,rhs;verbose=false,log=false);

@info "Computing setup pa 1x1 FMM $(M)DOFs freq$(Int(freq))"
LGM = LossyGlobalOuter(mesh,freq;fmm_on=true,depth=1,n=3,progress=false);
rhs = LGM.Ga*gmres(LGM.inner,(LGM.Dr*v0 - LGM.Nd'*gmres(LGM.Gv,LGM.Hv*v0;verbose=false));verbose=false);
paFMM = gmres(LGM,rhs;verbose=false,log=false);

xyzb = mesh.sources;
# Generating analytical solution
ρ,c,kₚ,kₐ,kₕ,kᵥ,τₐ,τₕ,ϕₐ,ϕₕ,η,μ = visco_thermal_constants(;freq=freq,S=1);
radius = 1.0;                                     # Radius of sphere_1m       [m]
coordinates = [radius*ones(M,1) acos.(xyzb[3,:]/radius)];
pasAN, v_rAN, v_thetaAN, v_rAN_A, v_thetaAN_A, v_rAN_V, v_thetaAN_V =
                BoundaryIntegralEquations.sphere_first_order(kₐ,c,ρ,radius,u₀,coordinates;S=1,kv=kᵥ);
ang_axis = coordinates[:,2]*180.0/pi;
perm = sortperm(ang_axis);

# Plotting pressure
K = 1;
gr(size=(600,500));
scatter(ang_axis[1:K:end],real.(pa[1:K:end]),label="BEM - 1n Iter",marker=:cross,markersize=2,color=:black,dpi=400);
scatter!(ang_axis[1:K:end],real.(paFMM[1:K:end]),label="BEM - 1n FMM",marker=:x,markersize=2,color=:red,dpi=400);
plot!(ang_axis[perm],real.(pasAN[perm]),label="Analytical",linewidth=1,color=:blue);
ylabel!("Re(Pa)");
title!("Frequency = $(freq) Hz");
xlabel!("Angle [deg]");
savefig("paGlobalFMM1x1_$(M)DOFs_$(Int(freq))Hz.png")


# Saving data
jldsave("resultsFMM1x1_$(M)DOFs_$(Int(freq))Hz.JLD2",
paFMM=paFMM,
pa=pa,
pasAN=pasAN,
ang_axis=ang_axis,
perm=perm,
);