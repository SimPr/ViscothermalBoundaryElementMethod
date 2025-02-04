# Save plots and solution for pa (1n-iter) and vv(1n-iter) + analytical reference
# Call with julia --project=.  ./examples/sphere_global_losses_1x1_resultsIGA.jl --freq 500 --compute_full_solution true --mesh_file "sphere_1m_coarser"

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
radius = 1.0;

#mesh = BoundaryIntegralEquations.load3dTriangularComsolMesh_SphereCorr_Disc(tri_mesh_file,radius;geometry_order=geometry_orders[2],
#                                                    physics_order=tri_physics_orders[2]);
mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[2],
                                                       physics_order=tri_physics_orders[5])

#= # Export geometry mesh and nodal normals/tangentials to paraview format (https://juliahub.com/ui/Packages/WriteVTK/D2v2J/1.12.0)
using WriteVTK
M  = size(mesh.sources,2);
numElem = size(mesh.physics_topology,2)
cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_TRIANGLE, mesh.physics_topology[:,i]) for i=1:numElem] # array that contains all the cells of the mesh 
vtkfile = vtk_grid("discMesh_$(M)DOFs", mesh.sources, cells)
vtkfile["n", VTKPointData()] = mesh.normals
vtkfile["t", VTKPointData()] = mesh.tangents
vtkfile["s", VTKPointData()] = mesh.sangents
outfiles = vtk_save(vtkfile)

using WriteVTK
M  = size(mesh.sources,2);
numElem = size(mesh.topology,2)
cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_TRIANGLE, mesh.topology[:,i]) for i=1:numElem] # array that contains all the cells of the mesh 
vtkfile = vtk_grid("contMesh_$(M)DOFs", mesh.coordinates, cells)
outfiles = vtk_save(vtkfile) =#

    

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

# iterative
@info "Assembling and solving iterative system"
pa = gmres(LGM,rhs;verbose=false);

# Generating analytical solution
ρ,c,kₚ,kₐ,kₕ,kᵥ,τₐ,τₕ,ϕₐ,ϕₕ,η,μ = visco_thermal_constants(;freq=float(freq),S=1);
radius = 1.0;                                     # Radius of sphere_1m       [m]
# polar and azimuth angle of spherical coordinates
theta = acos.(xyzb[3,:]./radius);
phi =  acos.(xyzb[1,:]./sqrt.(xyzb[1,:].^2 .+ xyzb[2,:].^2)).*sign.(xyzb[2,:]);
pasAN, v_rAN_V, v_thetaAN_V = BoundaryIntegralEquations.sphere_first_orderNew(u₀,kₚ,kᵥ,radius,theta,ρ,c);
ang_axis = theta*180.0/pi;
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
    v_r = LGM.Nd'*v;
    v_tang = v + LGM.Nd*v_r; # Computing the tangential velocity by substracting the normal information

    # recompute v_theta according to coordinate transformation from cartesian x,y,z to spherical theta
    v_theta = cos.(theta).*cos.(phi).*v_tang[0M+1:1M] + cos.(theta).*sin.(phi).*v_tang[1M+1:2M] - sin.(theta).*v_tang[2M+1:3M];
    #v_theta = sqrt.(v_tang[0M+1:1M].^2 + v_tang[1M+1:2M].^2 + v_tang[2M+1:3M].^2);
    #===========================================================================================
                                    Plotting solutions
    ===========================================================================================#
    # Plotting
    plt1 = scatter(ang_axis,real.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
    ylabel!("Re(Pa)"); plot!(ang_axis[perm],real.(pasAN[perm]),label="Analytical",linewidth=2,color=:blue);
    title!("Frequency = $(freq) Hz");
    plt2 = scatter(ang_axis,real.(v_r),label="BEM",marker=:cross,markersize=2,color=:black);
    ylabel!("Re(V_r)"); plot!(ang_axis[perm],real.(v_rAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
    plt3 = scatter(ang_axis,real.(v_theta),label="BEM",marker=:cross,markersize=2,color=:black);
    plot!(ang_axis[perm],real.(v_thetaAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
    xlabel!("Angle [deg]"); ylabel!("Re(V_theta)");
    plt4 = plot(plt1,plt2,plt3,layout=(3,1))
    savefig("allGlobal1x1b_disc_Real_$(M)DOFs_$(Int(freq))Hz.png")

    plt1 = scatter(ang_axis,imag.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
    ylabel!("Imag(Pa)"); plot!(ang_axis[perm],imag.(pasAN[perm]),label="Analytical",linewidth=2,color=:blue);
    title!("Frequency = $(freq) Hz");
    plt2 = scatter(ang_axis,imag.(v_r),label="BEM",marker=:cross,markersize=2,color=:black);
    ylabel!("Imag(V_r)"); plot!(ang_axis[perm],imag.(v_rAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
    plt3 = scatter(ang_axis,imag.(v_theta),label="BEM",marker=:cross,markersize=2,color=:black);
    plot!(ang_axis[perm],imag.(v_thetaAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
    xlabel!("Angle [deg]"); ylabel!("Imag(V_theta)");
    plt4 = plot(plt1,plt2,plt3,layout=(3,1))
    savefig("allGlobal1x1b_disc_Imag_$(M)DOFs_$(Int(freq))Hz.png")

    plt1 = scatter(ang_axis,abs.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
    ylabel!("|Pa|"); plot!(ang_axis[perm],abs.(pasAN[perm]),label="Analytical",linewidth=2,color=:blue);
    title!("Frequency = $(freq) Hz");
    plt2 = scatter(ang_axis,abs.(v_r),label="BEM",marker=:cross,markersize=2,color=:black);
    ylabel!("|V_r|"); plot!(ang_axis[perm],abs.(v_rAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
    plt3 = scatter(ang_axis,abs.(v_theta),label="BEM",marker=:cross,markersize=2,color=:black);
    plot!(ang_axis[perm],abs.(v_thetaAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
    xlabel!("Angle [deg]"); ylabel!("|V_theta|");
    plt4 = plot(plt1,plt2,plt3,layout=(3,1))
    savefig("allGlobal1x1b_disc_Abs_$(M)DOFs_$(Int(freq))Hz.png")

    # Saving data
    jldsave("results1x1b_disc_$(M)DOFs_$(Int(freq))Hz.JLD2", 
    pa=pa,
    pasAN=pasAN,
    v_rAN_V=v_rAN_V,
    v_thetaAN_V=v_thetaAN_V,
    ang_axis=ang_axis,
    v_r=v_r,
    v_theta=v_theta,
    hist_dpa=hist_dpa,
    )

else
    # Saving data
    jldsave("results1x1b_disc_$(M)DOFs_$(Int(freq))Hz.JLD2", 
    pa=pa,
    pasAN=pasAN,
    ang_axis=ang_axis,
    perm=perm,
    );
end

function errorCalc(calc,ref,M) # calculate relative root mean square error

    eps = (sum(abs.(calc-ref).^2)/M).^(1/2)
    refNorm = (sum(abs.(ref).^2)/M).^(1/2);
    epsRel = eps/refNorm;

    return epsRel
end

#v_abs = sqrt.(v_r.^2 + v_theta.^2);
v_glob_abs = sqrt.(v[0M+1:1M].^2 + v[1M+1:2M].^2 + v[2M+1:3M].^2);
v_absAN = sqrt.(v_rAN_V.^2 + v_thetaAN_V.^2);

eps_pa = errorCalc(pa,pasAN,M)
eps_vvr = errorCalc(v_r,v_rAN_V,M)
eps_vvtheta = errorCalc(v_theta,v_thetaAN_V,M)
eps_vv = errorCalc(v_glob_abs,v_absAN,M)

eps_paDisc = errorCalc(pa,pasAN,M)
eps_vvrDisc = errorCalc(v_r,v_rAN_V,M)
eps_vvthetaDisc = errorCalc(v_theta,v_thetaAN_V,M)
eps_vvDisc = errorCalc(v_glob_abs,v_absAN,M)

using DelimitedFiles
results = zeros(Float64, M,6); 
results[:,1] = real.(pa)
results[:,2] = imag.(pa)
results[:,4] = real.(v_r)
results[:,5] = imag.(v_r)
results[:,7] = real.(v_theta)
results[:,8] = imag.(v_theta)

#writedlm( "results_wing.csv",  results, ',')
writedlm( "results_disc_sphere_$(M)DOFs_$(Int(freq))Hz.txt",  results, '\t')
writedlm( "coordinates_disc_sphere_$(M)DOFs.txt",  mesh.sources, '\t')
writedlm( "elements_disc_sphere_$(M)DOFs.txt",  mesh.physics_topology, '\t')


#==========================================================================================
                            Export corrected mesh
==========================================================================================#
mesh_path = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes");
mesh_files =  ["sphere_1m_coarser","sphere_1m_coarse","sphere_1m","sphere_1m_fine","sphere_1m_finer","sphere_1m_4p5k","sphere_1m_extremely_fine","sphere_1m_35k","sphere_1m_77k"];
mesh_file = mesh_files[8]
tri_mesh_file = joinpath(mesh_path,mesh_file);
radius = 1.0;
mesh = BoundaryIntegralEquations.load3dTriangularComsolMesh_SphereCorr(tri_mesh_file,radius;geometry_order=geometry_orders[2],
                                                    physics_order=tri_physics_orders[2]);

using WriteVTKmesh
M  = size(mesh.sources,2);
numElem = size(mesh.topology,2)
cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_TRIANGLE, mesh.topology[:,i]) for i=1:numElem] # array that contains all the cells of the mesh 
vtkfile = vtk_grid("$(mesh_file)_corr", mesh.coordinates, cells)
outfiles = vtk_save(vtkfile)