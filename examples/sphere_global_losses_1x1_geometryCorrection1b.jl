# Save plots and solution for paDense (1n-iter) and vv(1n-iter) + analytical reference
# Call with julia --project=.  ./examples/sphere_global_losses_1x1_results1b.jl --freq 1000 --compute_full_solution true --mesh_file "sphere_1m_fine"

using ArgParse
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--freq"
            help = "frequency"
            default = 1000.0
            arg_type = Float64
        "--compute_full_solution"
            help = raw"computing full solution"
            default = true
            arg_type = Bool
        "--mesh_file"
            help = "mesh file"
            default = "sphere_1m_fine"
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

# iter
@info "Solving iterative system"
pa = gmres(LGM,rhs;verbose=false);

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
    plt4 = plot(plt1,plt2,plt3,layout=(3,1))
    savefig("allGlobal1x1b_Real_$(M)DOFs_$(Int(freq))Hz.png")

    plt1 = scatter(ang_axis,imag.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
    ylabel!("Imag(Pa)"); plot!(ang_axis[perm],imag.(pasAN[perm]),label="Analytical",linewidth=2,color=:blue);
    title!("Frequency = $(freq) Hz");
    plt2 = scatter(ang_axis,imag.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black);
    ylabel!("Imag(Vn)"); plot!(ang_axis[perm],imag.(v_rAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
    plt3 = scatter(ang_axis,imag.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black);
    plot!(ang_axis[perm],imag.(v_thetaAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
    xlabel!("Angle [deg]"); ylabel!("Imag(Vt)");
    plt4 = plot(plt1,plt2,plt3,layout=(3,1))
    savefig("allGlobal1x1b_Imag_$(M)DOFs_$(Int(freq))Hz.png")

    plt1 = scatter(ang_axis,abs.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
    ylabel!("|Pa|"); plot!(ang_axis[perm],abs.(pasAN[perm]),label="Analytical",linewidth=2,color=:blue);
    title!("Frequency = $(freq) Hz");
    plt2 = scatter(ang_axis,abs.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black);
    ylabel!("|Vn|"); plot!(ang_axis[perm],abs.(v_rAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
    plt3 = scatter(ang_axis,abs.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black);
    plot!(ang_axis[perm],abs.(v_thetaAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
    xlabel!("Angle [deg]"); ylabel!("|Vt|");
    plt4 = plot(plt1,plt2,plt3,layout=(3,1))
    savefig("allGlobal1x1b_Abs_$(M)DOFs_$(Int(freq))Hz.png")

    # Saving data
    jldsave("results1x1b_$(M)DOFs_$(Int(freq))Hz.JLD2", 
    pa=pa,
    pasAN=pasAN,
    v_rAN_V=v_rAN_V,
    v_thetaAN_V=v_thetaAN_V,
    ang_axis=ang_axis,
    v_n0=v_n0,
    vt_sum=vt_sum,
    hist_dpa=hist_dpa,
    )

else
    # Saving data
    jldsave("results1x1b_$(M)DOFs_$(Int(freq))Hz.JLD2", 
    pa=pa,
    pasAN=pasAN,
    ang_axis=ang_axis,
    perm=perm,
    );
end

#===========================================================================================
                        Alternative analytical solution and computation of vt_sum
===========================================================================================#
function sph_hankel(z)
    # returns spherical hankel function 1st kind 1st order and its derivative
    h1  = -exp.(1im*z)./z.^2 .*(z .+ 1im)
    dh1 = exp.(1im*z)./z.^3 .*(2.0.*z .+ 1im.*(-z.^2 .+2.0))
    return h1,dh1
end

function compute_a1(v0,k,kv,a)
    # returns coefficient a1 needed for analytical solution
    a1=-v0*k^2*a*exp(-1im*k*a)*(3*kv*a+3*1im-1im*kv^2*a^2)/(3*1im*(kv^2*(k^2*a^2-2.0)-k^2+1im*kv*k*a*(k+2*kv)))
    return a1
end

function sphere_first_order2(v0,k,kv,a,theta,rho,c)
    # returns analytical solution for sphere transversely oscillating in viscous fluid
    a1 = compute_a1(v0,k,kv,a)
    h1_ka,dh1_ka = sph_hankel(k*a)
    pa = -3*rho*k*c*a1*h1_ka.*cos.(theta)
    vv_r = (v0 - 3*1im*k*a1*dh1_ka).*cos.(theta)
    vv_theta = (-v0 + 3*1im/a*a1*h1_ka).*sin.(theta)
    return pa, vv_r, vv_theta
end

pasAN2, v_rAN_V2, v_thetaAN_V2 = sphere_first_order2(u₀,kₚ,kᵥ,radius,coordinates[:,2],ρ,c)
# polar and azimuth angle of spherical coordinates
theta = acos.(xyzb[3,:]./radius)
phi =  acos.(xyzb[1,:]./sqrt.(xyzb[1,:].^2 .+ xyzb[2,:].^2)).*xyzb[2,:]./abs.(xyzb[2,:])

# recompute vt_sum according to coordinate transformation from cartesian x,y,z to spherical theta
vt_sum = cos.(theta).*cos.(phi).*v_t[0M+1:1M] + cos.(theta).*sin.(phi).*v_t[1M+1:2M] - sin.(theta).*v_t[2M+1:3M]

plt1 = scatter(ang_axis,real.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
ylabel!("Re(Pa)"); plot!(ang_axis[perm],real.(pasAN2[perm]),label="Analytical2",linewidth=2,color=:blue);
title!("Frequency = $(freq) Hz");
plt2 = scatter(ang_axis,real.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black);
ylabel!("Re(Vn)"); plot!(ang_axis[perm],real.(v_rAN_V2[perm]),label="Analytical",linewidth=2,color=:blue);
plt3 = scatter(ang_axis,real.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black);
plot!(ang_axis[perm],real.(v_thetaAN_V2[perm]),label="Analytical2",linewidth=2,color=:blue);
xlabel!("Angle [deg]"); ylabel!("Re(Vt)");
plt4 = plot(plt1,plt2,plt3,layout=(3,1))
savefig("allGlobal1x1b_Real_$(M)DOFs_$(Int(freq))Hz2.png")

plt1 = scatter(ang_axis,imag.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
ylabel!("Imag(Pa)"); plot!(ang_axis[perm],imag.(pasAN2[perm]),label="Analytical2",linewidth=2,color=:blue);
title!("Frequency = $(freq) Hz");
plt2 = scatter(ang_axis,imag.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black);
ylabel!("Imag(Vn)"); plot!(ang_axis[perm],imag.(v_rAN_V2[perm]),label="Analytical2",linewidth=2,color=:blue);
plt3 = scatter(ang_axis,imag.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black);
plot!(ang_axis[perm],imag.(v_thetaAN_V2[perm]),label="Analytical2",linewidth=2,color=:blue);
xlabel!("Angle [deg]"); ylabel!("Imag(Vt)");
plt4 = plot(plt1,plt2,plt3,layout=(3,1))
savefig("allGlobal1x1b_Imag_$(M)DOFs_$(Int(freq))Hz2.png")

plt1 = scatter(ang_axis,abs.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
ylabel!("|Pa|"); plot!(ang_axis[perm],abs.(pasAN2[perm]),label="Analytical2",linewidth=2,color=:blue);
title!("Frequency = $(freq) Hz");
plt2 = scatter(ang_axis,abs.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black);
ylabel!("|Vn|"); plot!(ang_axis[perm],abs.(v_rAN_V2[perm]),label="Analytical2",linewidth=2,color=:blue);
plt3 = scatter(ang_axis,abs.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black);
plot!(ang_axis[perm],abs.(v_thetaAN_V2[perm]),label="Analytical2",linewidth=2,color=:blue);
xlabel!("Angle [deg]"); ylabel!("|Vt|");
plt4 = plot(plt1,plt2,plt3,layout=(3,1))
savefig("allGlobal1x1b_Abs_$(M)DOFs_$(Int(freq))Hz2.png")

jldsave("results1x1b_$(M)DOFs_$(Int(freq))Hz2.JLD2", 
       pa=pa,
       pasAN=pasAN2,
       v_rAN_V=v_rAN_V2,
       v_thetaAN_V=v_thetaAN_V2,
       ang_axis=ang_axis,
       v_n0=v_n0,
       vt_sum=-vt_sum,
       hist_dpa=hist_dpa,
       )


#==========================================================================================
                            Compute gometry error and redo Simulation for corrected mesh
==========================================================================================#

r = sqrt.(xyzb[1,:].^2+xyzb[2,:].^2+xyzb[3,:].^2);

function errorCalc(calc,ref,M)

    eps = (sum(abs.(calc.-ref).^2)/M).^(1/2)
    refNorm = (sum(abs.(ref).^2)/M).^(1/2);
    epsRel = eps/refNorm;
    return epsRel
end

epsRel = errorCalc(r,radius,2);
scatter(ang_axis[perm],(radius.-r)[perm],label="p2 $(M)DOFs epsRel=$(round(epsRel; digits = 7))",marker=:cross,markersize=2,color=:black,dpi=400); # for quadratic physics discretization
title!("Geometry error");
xlabel!("Angle [deg]");
ylabel!("r_ref - r_num [m]");
savefig("geometryError.png")


pa_error = errorCalc(pa,pasAN,M) # discretization error in surface nodes compared to analytical solution
vn_error = errorCalc(v_n0,v_rAN_V,M)
vt_error = errorCalc(vt_sum,v_thetaAN_V,M)


###############################


meshCorr = BoundaryIntegralEquations.load3dTriangularComsolMesh_SphereCorr(tri_mesh_file,radius;geometry_order=geometry_orders[2],
                                                    physics_order=tri_physics_orders[2]);

#==========================================================================================
                            Creating excitation vector
==========================================================================================#
xyzbCorr = meshCorr.sources;
#===========================================================================================
                        BEM matrix assembly and (iterative) solution of the 1-variable system
===========================================================================================#

LGM = LossyGlobalOuter(meshCorr,freq;fmm_on=false,depth=1,n=3,progress=false);

@info "Computing RHS"
rhs = LGM.Ga*gmres(LGM.inner,(LGM.Dr*v0 - LGM.Nd'*gmres(LGM.Gv,LGM.Hv*v0;verbose=false));verbose=false);

# iter
@info "Solving iterative system"
pa = gmres(LGM,rhs;verbose=false);

# Generating analytical solution
coordinates = [radius*ones(M,1) acos.(xyzbCorr[3,:]/radius)];
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
    plt4 = plot(plt1,plt2,plt3,layout=(3,1))
    savefig("allGlobal1x1bCorr_Real_$(M)DOFs_$(Int(freq))Hz.png")

    plt1 = scatter(ang_axis,imag.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
    ylabel!("Imag(Pa)"); plot!(ang_axis[perm],imag.(pasAN[perm]),label="Analytical",linewidth=2,color=:blue);
    title!("Frequency = $(freq) Hz");
    plt2 = scatter(ang_axis,imag.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black);
    ylabel!("Imag(Vn)"); plot!(ang_axis[perm],imag.(v_rAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
    plt3 = scatter(ang_axis,imag.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black);
    plot!(ang_axis[perm],imag.(v_thetaAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
    xlabel!("Angle [deg]"); ylabel!("Imag(Vt)");
    plt4 = plot(plt1,plt2,plt3,layout=(3,1))
    savefig("allGlobal1x1bCorr_Imag_$(M)DOFs_$(Int(freq))Hz.png")

    plt1 = scatter(ang_axis,abs.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
    ylabel!("|Pa|"); plot!(ang_axis[perm],abs.(pasAN[perm]),label="Analytical",linewidth=2,color=:blue);
    title!("Frequency = $(freq) Hz");
    plt2 = scatter(ang_axis,abs.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black);
    ylabel!("|Vn|"); plot!(ang_axis[perm],abs.(v_rAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
    plt3 = scatter(ang_axis,abs.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black);
    plot!(ang_axis[perm],abs.(v_thetaAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
    xlabel!("Angle [deg]"); ylabel!("|Vt|");
    plt4 = plot(plt1,plt2,plt3,layout=(3,1))
    savefig("allGlobal1x1bCorr_Abs_$(M)DOFs_$(Int(freq))Hz.png")

    # Saving data
    jldsave("results1x1bCorr_$(M)DOFs_$(Int(freq))Hz.JLD2", 
    pa=pa,
    pasAN=pasAN,
    v_rAN_V=v_rAN_V,
    v_thetaAN_V=v_thetaAN_V,
    ang_axis=ang_axis,
    v_n0=v_n0,
    vt_sum=vt_sum,
    hist_dpa=hist_dpa,
    )

else
    # Saving data
    jldsave("results1x1bCorr_$(M)DOFs_$(Int(freq))Hz.JLD2", 
    pa=pa,
    pasAN=pasAN,
    ang_axis=ang_axis,
    perm=perm,
    );
end

# Geometry error
rCorr = sqrt.(xyzbCorr[1,:].^2+xyzbCorr[2,:].^2+xyzbCorr[3,:].^2);

epsRelCorr = errorCalc(rCorr,radius,2);
scatter(ang_axis[perm],(radius.-rCorr)[perm],label="p2 $(M)DOFs epsRelCorr=$(round(epsRelCorr; digits = 7))",marker=:cross,markersize=2,color=:black,dpi=400); # for quadratic physics discretization
title!("Geometry error");
xlabel!("Angle [deg]");
ylabel!("r_ref - r_numCorr [m]");
savefig("geometryErrorCorr.png")


pa_errorCorr = errorCalc(pa,pasAN,M) # discretization error in surface nodes compared to analytical solution
vn_errorCorr = errorCalc(v_n0,v_rAN_V,M)
vt_errorCorr = errorCalc(vt_sum,v_thetaAN_V,M)


#######################################

using DelimitedFiles

data_file = "results1x1bCorr_$(2642)DOFs_$(Int(1000))Hz.JLD2"
f = jldopen(data_file)

pa = f["pa"]
v_n0 = f["v_n0"]
vt_sum = f["vt_sum"]
ang_axis = f["ang_axis"]

perm = sortperm(ang_axis)
writedlm( "ang_axisCorr.csv",  ang_axis, ',')
writedlm( "ang_axisPermCorr.csv",  ang_axis[perm], ',')
writedlm( "paRealCorr.csv",  real.(pa), ',')
writedlm( "paImagCorr.csv", imag.(pa), ',')
writedlm( "paAbsCorr.csv",  abs.(pa), ',')
writedlm( "vnRealCorr.csv",  real.(v_n0), ',')
writedlm( "vnImagCorr.csv",  imag.(v_n0), ',')
writedlm( "vnAbsCorr.csv",  abs.(v_n0), ',')
writedlm( "vtRealCorr.csv",  real.(vt_sum), ',')
writedlm( "vtImagCorr.csv",  imag.(vt_sum), ',')
writedlm( "vtAbsCorr.csv",  abs.(vt_sum), ',')
close(f)


