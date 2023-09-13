# Better name would be sphere_global_losses_1x1_results1
# Save plots and solution for paDense (1n-direct) and vv(1n-direct) + analytical reference
# Call with julia --project=.  ./examples/sphere_global_losses_1x1_runTimeComp1.jl --freq 500 --compute_full_solution true --mesh_file "sphere_1m_coarser"

using ArgParse
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--freq"
            help = "frequency"
            default = 100.0
            arg_type = Float64
        "--compute_full_solution"
            help = raw"computing full solution"
            default = false
            arg_type = Bool
        "--mesh_file"
            help = "mesh file"
            default = "sphere_1m_coarser"
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
compute_full_solution = intputArguments["compute_full_solution"]
mesh_file = intputArguments["mesh_file"]

#==========================================================================================
                            Adding Related Packages
==========================================================================================#
using LinearAlgebra
using BoundaryIntegralEquations
using Plots
using IterativeSolvers
using JLD2


#=============================ß============================================================
                                Loading Mesh
==========================================================================================#
geometry_orders     = [:linear,:quadratic];
tri_physics_orders  = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic];
#bool_reconstruction = true;
# Triangular Meshes

#mesh_files =  ["sphere_1m_coarser"]#,"sphere_1m_coarse","sphere_1m","sphere_1m_fine","sphere_1m_finer","sphere_1m_4p5k","sphere_1m_extremely_fine","sphere_1m_35k","sphere_1m_77k"];
mesh_path = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes");

#mesh_file = mesh_files[1]
tri_mesh_file = joinpath(mesh_path,mesh_file);
mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[2],
                                                    physics_order=tri_physics_orders[2]);
#==========================================================================================
                                Setting up constants
==========================================================================================#
#freq   = 1000.0;                                   # Frequency                 [Hz]
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

LGM = LossyGlobalOuter(mesh,freq;fmm_on=false,depth=1,n=3,progress=false);

@info "Computing RHS"
rhs = LGM.Ga*gmres(LGM.inner,(LGM.Dr*v0 - LGM.Nd'*gmres(LGM.Gv,LGM.Hv*v0;verbose=false));verbose=false);

# dense
@info "Assembling and solving dense system"
LGM_dense  = BoundaryIntegralEquations._full1_new(LGM);
pa = LGM_dense\rhs;

# Generating analytical solution
ρ,c,kₚ,kₐ,kₕ,kᵥ,τₐ,τₕ,ϕₐ,ϕₕ,η,μ = visco_thermal_constants(;freq=freq,S=1);
radius = 1.0;                                     # Radius of sphere_1m       [m]
coordinates = [radius*ones(M,1) acos.(xyzb[3,:]/radius)];
pasAN, v_rAN, v_thetaAN, v_rAN_A, v_thetaAN_A, v_rAN_V, v_thetaAN_V =
                BoundaryIntegralEquations.sphere_first_order(kₐ,c,ρ,radius,u₀,coordinates;S=1,kv=kᵥ);
ang_axis = coordinates[:,2]*180.0/pi;
perm = sortperm(ang_axis);

if compute_full_solution == true
    #===========================================================================================
                                    Reconstructing unknowns
    ===========================================================================================#
    @info "Reconstructing unknowns"
    tmp1 =  gmres(LGM.Gh,LGM.Hh*pa;verbose=false);
    dpa,hist_dpa  = gmres(LGM.Ga,LGM.Ha*pa;verbose=false,log=true); # <- This is the bottleneck...
    v  = v0 - (LGM.mu_a*LGM.Dc*pa + LGM.mu_h*LGM.Nd*tmp1 - LGM.phi_a*LGM.Nd*dpa);
    # Local components of the viscous velocity on the boundary
    v_n0 = LGM.Nd'*v;
    v_t = v + LGM.Nd*v_n0; # Computing the tangential velocity by substracting the normal information
    vt_sum = sqrt.(v_t[0M+1:1M].^2 + v_t[1M+1:2M].^2 + v_t[2M+1:3M].^2);
    #===========================================================================================
                                    Plotting solutions
    ===========================================================================================#
    # Plotting
    plt1 = scatter(ang_axis,real.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
    ylabel!("Re(Pa)"); plot!(ang_axis[perm],real.(pasAN[perm]),label="Analytical",linewidth=2,color=:blue);
    title!("Frequency = $(freq) Hz");
    plt2 = scatter(ang_axis,real.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black);
    ylabel!("Re(Vn)"); plot!(ang_axis[perm],real.(v_rAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
    plt3 = scatter(ang_axis,real.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black);
    plot!(ang_axis[perm],real.(v_thetaAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
    xlabel!("Angle [deg]"); ylabel!("Re(Vt)");
    plt4 = plot(plt1,plt2,plt3,layout=(3,1));
    savefig("allGlobal1x1_Real_$(M)DOFs_$(Int(freq))Hz.png")

    plt1 = scatter(ang_axis,imag.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
    ylabel!("Imag(Pa)"); plot!(ang_axis[perm],imag.(pasAN[perm]),label="Analytical",linewidth=2,color=:blue);
    title!("Frequency = $(freq) Hz");
    plt2 = scatter(ang_axis,imag.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black);
    ylabel!("Imag(Vn)"); plot!(ang_axis[perm],imag.(v_rAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
    plt3 = scatter(ang_axis,imag.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black);
    plot!(ang_axis[perm],imag.(v_thetaAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
    xlabel!("Angle [deg]"); ylabel!("Imag(Vt)");
    plt4 = plot(plt1,plt2,plt3,layout=(3,1));
    savefig("allGlobal1x1_Imag_$(M)DOFs_$(Int(freq))Hz.png")

    plt1 = scatter(ang_axis,abs.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
    ylabel!("|Pa|"); plot!(ang_axis[perm],abs.(pasAN[perm]),label="Analytical",linewidth=2,color=:blue);
    title!("Frequency = $(freq) Hz");
    plt2 = scatter(ang_axis,abs.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black);
    ylabel!("|Vn|"); plot!(ang_axis[perm],abs.(v_rAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
    plt3 = scatter(ang_axis,abs.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black);
    plot!(ang_axis[perm],abs.(v_thetaAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
    xlabel!("Angle [deg]"); ylabel!("|Vt|");
    plt4 = plot(plt1,plt2,plt3,layout=(3,1));
    savefig("allGlobal1x1_Abs_$(M)DOFs_$(Int(freq))Hz.png")

    # Saving data
    jldsave("results1x1_$(M)DOFs_$(Int(freq))Hz.JLD2", 
    paDense=pa,
    pasAN=pasAN,
    v_rAN_V=v_rAN_V,
    v_thetaAN_V=v_thetaAN_V,
    ang_axis=ang_axis,
    v_n0Dense=v_n0,
    vt_sumDense=vt_sum,
    hist_dpaDense=hist_dpa,
    )

else
    # Saving data
    jldsave("results1x1_$(M)DOFs_$(Int(freq))Hz.JLD2", 
    paDense=pa,
    pasAN=pasAN,
    ang_axis=ang_axis,
    perm=perm,
    );
end