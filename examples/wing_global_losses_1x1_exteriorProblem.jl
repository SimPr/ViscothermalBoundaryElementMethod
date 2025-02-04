# Save plots and solution for pa (1n-iter) and vv(1n-iter) + analytical reference
# Call with julia --project=.  ./examples/wing_global_losses_1x1_resultsIGA.jl --freq 500 --compute_full_solution true --mesh_file "sphere_1m_coarser"

using ArgParse
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--freq"
            help = "frequency"
            default = 500.0
            arg_type = Float64
        "--compute_full_solution"
            help = raw"computing full solution"
            default = true
            arg_type = Bool
        "--mesh_file"
            help = "mesh file"
            default = "wing"
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
# Triangular Meshes

#mesh_files =  ["sphere_1m_coarser"]#,"sphere_1m_coarse","sphere_1m","sphere_1m_fine","sphere_1m_finer","sphere_1m_4p5k","sphere_1m_extremely_fine","sphere_1m_35k","sphere_1m_77k"];
mesh_path = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes");

#mesh_file = mesh_files[1]
tri_mesh_file = joinpath(mesh_path,mesh_file);
mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[2],
                                                    physics_order=tri_physics_orders[5],int_ext=-1);
#==========================================================================================
                                Setting up constants
==========================================================================================#
#freq   = 1000.0;                                   # Frequency                 [Hz]
#==========================================================================================
                            Creating excitation vector
==========================================================================================#
xyzb = mesh.sources;
M  = size(xyzb,2);
u₀ = 0.01; # [m/s]
v0y = zeros(M);

#applying the bc when using the comsol entities
bc_entities = [2,3].-1
bc_element_id = Bool.(sum(bc_entities' .∈ mesh.entities,dims=2))[:]
bc_topology = mesh.physics_topology[:,bc_element_id]
bc_node_id = sort(unique(bc_topology)) # gives the same as findall(x->x==true, isapprox.(xyzb[1,:], 0.; atol=1E-10))
#using DelimitedFiles
#writedlm( "wing_bc_element_id.txt",  bc_element_id, '\t')

v0y[bc_node_id] .= u₀;
v0 = [zeros(M); v0y; zeros(M)];

using WriteVTK
numElem = size(mesh.physics_topology,2)
cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_TRIANGLE, mesh.physics_topology[:,i]) for i=1:numElem] # array that contains all the cells of the mesh 
vtkfile = vtk_grid("wing_v0_$(M)physDOFs", mesh.sources, cells)
#numElem = size(bc_topology,2)
#numNode = size(bc_node_id,1)
#cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_TRIANGLE, bc_topology[:,i]) for i=1:numElem] # array that contains all the cells of the mesh 
#vtkfile = vtk_grid("wing_v0_$(M)physDOFs", mesh.sources[:,bc_node_id], cells)
v0_BC = reshape(v0, (M,3))'
vtkfile["v0", VTKPointData()] = v0_BC#[:,bc_node_id]
outfiles = vtk_save(vtkfile)

#= data_file = "results1x1b_$(M)DOFs_$(Int(freq))Hz.JLD2"
f = jldopen(data_file)
pa=f["pa"];
#hist_pa=f["hist_pa"];
v=f["v"];
#v_r=f["v_r"];
#v_theta=f["v_theta"];
#hist_dpa=f["hist_dpa"];
close(f)

v_abs = sqrt.(v[0M+1:1M].^2 + v[1M+1:2M].^2 + v[2M+1:3M].^2);
log_vv = log10.(abs.(v_abs));

LGM = LossyGlobalOuter(mesh,freq;fmm_on=false,depth=1,n=3,progress=false,int_ext=-1);
p_h = -LGM.tau_a/LGM.tau_h.*pa;
log_ph = log10.(abs.(p_h));

p0 = 2*10^(-5);
pa_rms = sqrt.(0.5*real.(pa.*conj.(pa)));
SPL = 20*log10.(pa_rms./p0);

using DelimitedFiles
results = zeros(Float64, M,6); 
results[:,1] = xyzb[1,:];
results[:,2] = xyzb[2,:];
results[:,3] = xyzb[3,:];
results[:,4] = SPL;
results[:,5] = log_ph;
results[:,6] = log_vv;

#writedlm( "results_wing.csv",  results, ',')
writedlm( "results_wing_extfine.txt",  results, '\t')

writedlm( "coordinates_wing_extfine.txt",  mesh.sources, '\t')
writedlm( "elements_wing_extfine.txt",  mesh.physics_topology, '\t')
writedlm( "SPL_wing_extfine.txt",  SPL, '\t')
writedlm( "log_ph_wing_extfine.txt",  log_ph, '\t')
writedlm( "log_vv_wing_extfine.txt",  log_vv, '\t')

SPL = readdlm("SPL_wing_extfine.txt", '\t', Float64, '\n')
log_ph = readdlm("log_ph_wing_extfine.txt", '\t', Float64, '\n')
log_vv = readdlm("log_vv_wing_extfine.txt", '\t', Float64, '\n')

# Export physics mesh and nodal normals/tangentials and nodal results to paraview format (https://juliahub.com/ui/Packages/WriteVTK/D2v2J/1.12.0)
using WriteVTK
numElem = size(mesh.physics_topology,2)
cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_TRIANGLE, mesh.physics_topology[:,i]) for i=1:numElem] # array that contains all the cells of the mesh 
vtkfile = vtk_grid("wing_$(M)physDOFs", mesh.sources, cells)
vtkfile["n", VTKPointData()] = mesh.normals
vtkfile["t", VTKPointData()] = mesh.tangents
vtkfile["s", VTKPointData()] = mesh.sangents
vtkfile["SPL", VTKPointData()] = SPL
vtkfile["log10ph", VTKPointData()] = log_ph
vtkfile["log10vv", VTKPointData()] = log_vv
outfiles = vtk_save(vtkfile) 
=#

using WriteVTK
cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_TRIANGLE, mesh.topology[:,i]) for i=1:size(mesh.topology,2)] # array that contains all the cells of the mesh 
vtkfile = vtk_grid("wing_$(size(mesh.sources,2))physDOFs", mesh.coordinates, cells)
outfiles = vtk_save(vtkfile) 
#===========================================================================================
                        BEM matrix assembly and (iterative) solution of the 1-variable system
===========================================================================================#

LGM = LossyGlobalOuter(mesh,freq;fmm_on=false,depth=1,n=3,progress=false,int_ext=-1);

@info "Computing RHS"
rhs = LGM.Ga*gmres(LGM.inner,(LGM.Dr*v0 - LGM.Nd'*gmres(LGM.Gv,LGM.Hv*v0;verbose=false));verbose=false);

# iterative
@info "Assembling and solving iterative system"
pa,hist_pa = gmres(LGM,rhs;verbose=true,log=true);

p0 = 2*10^(-5);
pa_rms = sqrt.(0.5*real.(pa.*conj.(pa)));
SPL = 20*log10.(pa_rms./p0);

ph = -LGM.tau_a/LGM.tau_h.*pa;
log_ph = log10.(abs.(ph));


if compute_full_solution == true
    #===========================================================================================
                                    Reconstructing unknowns
    ===========================================================================================#
    @info "Reconstructing unknowns"
    tmp1 =  gmres(LGM.Gh,LGM.Hh*pa;verbose=false);
    dpa,hist_dpa  = gmres(LGM.Ga,LGM.Ha*pa;verbose=true,log=true); # <- This is the bottleneck...
    v  = v0 - (LGM.mu_a*LGM.Dc*pa + LGM.mu_h*LGM.Nd*tmp1 - LGM.phi_a*LGM.Nd*dpa);
    # Local components of the viscous velocity on the boundary
    v_r = LGM.Nd'*v;
    v_tang = v + LGM.Nd*v_r; # Computing the tangential velocity by substracting the normal information

    # recompute v_theta according to coordinate transformation from cartesian x,y,z to spherical theta
    #v_theta = cos.(theta).*cos.(phi).*v_tang[0M+1:1M] + cos.(theta).*sin.(phi).*v_tang[1M+1:2M] - sin.(theta).*v_tang[2M+1:3M];
    v_theta = sqrt.(v_tang[0M+1:1M].^2 + v_tang[1M+1:2M].^2 + v_tang[2M+1:3M].^2);

    v_abs = sqrt.(v[0M+1:1M].^2 + v[1M+1:2M].^2 + v[2M+1:3M].^2);
    log_vv = log10.(abs.(v_abs));

    # Saving data
    jldsave("results1x1b_$(M)DOFs_$(Int(freq))Hz.JLD2", 
    pa=pa,
    SPL=SPL,
    ph=ph,
    log_ph=log_ph,
    hist_pa=hist_pa,
    v=v,
    log_vv=log_vv,
    v_r=v_r,
    v_theta=v_theta,
    hist_dpa=hist_dpa,
    coordinates=mesh.sources,
    elements=mesh.physics_topology,
    )
    using DelimitedFiles
    writedlm( "coordinates_wing.txt",  mesh.sources, '\t')
    writedlm( "elements_wing.txt",  mesh.physics_topology, '\t')
    writedlm( "SPL_wing.txt",  SPL, '\t')
    writedlm( "log_ph_wing.txt",  log_ph, '\t')
    writedlm( "log_vv_wing.txt",  log_vv, '\t')

else
    # Saving data
    jldsave("results1x1b_$(M)DOFs_$(Int(freq))Hz.JLD2", 
    pa=pa,
    SPL=SPL,
    ph=ph,
    log_ph=log_ph,
    hist_pa=hist_pa,
    coordinates=mesh.sources,
    elements=mesh.physics_topology,
    );
end

#===========================================================================================
                                Plotting solutions
===========================================================================================#

#= using MeshViz
import WGLMakie as wgl
# simple_mesh = create_simple_mesh(mesh);
wgl.set_theme!(resolution=(800, 800));
tri_bc_ents = [2,3].-1 # top abd bottom surface of wing
# Plotting boundary conditions
tri_simple_mesh = create_bc_simple_mesh(mesh,tri_bc_ents,false);
viz(tri_simple_mesh;showfacets=true,alpha=0.5);
tri_simple_bc   = create_bc_simple_mesh(mesh,tri_bc_ents[1]);
viz!(tri_simple_bc;showfacets=true,color=:red);
tri_simple_bc   = create_bc_simple_mesh(mesh,tri_bc_ents[2]);
viz!(tri_simple_bc;showfacets=true,color=:blue) 
tri_simple_points   = BoundaryIntegralEquations.create_points(mesh,mesh.sources[:,bc_node_id])
viz!(tri_simple_points;pointsize=10,colorscheme=:jet1)    

# plotting SPL, pressure, and BC
p0 = 2*10^(-5);
pa_rms = sqrt.(0.5*real.(pa.*conj.(pa)));
SPL = 20*log10.(pa_rms./p0);
wgl.set_theme!(resolution=(1000, 1000));
data_mesh,data_viz = create_vizualization_data(mesh,SPL);
fig, ax, hm = viz(data_mesh;showfacets=false,color=data_viz,colorscheme=:jet1)
wgl.Colorbar(fig[1,2], colorrange = (minimum(data_viz),maximum(data_viz)), colormap = :jet1, label="SPL (dB)")
wgl.save("SPL1x1b_$(M)DOFs_$(Int(freq))Hz.html", fig)

data_mesh,data_viz = create_vizualization_data(mesh,abs.(pa));
fig2, ax, hm = viz(data_mesh;showfacets=false, color=data_viz,colorscheme=:jet1)
wgl.Colorbar(fig2[1,2], colorrange = (minimum(data_viz),maximum(data_viz)), colormap = :jet1, label="|p_a| (Pa)");
wgl.save("pa_abs_1x1b_$(M)DOFs_$(Int(freq))Hz.html", fig2)

v0_abs = sqrt.(v0[0M+1:1M].^2 + v0[1M+1:2M].^2 + v0[2M+1:3M].^2); 
data_mesh,data_viz = create_vizualization_data(mesh,v0_abs);
fig3, ax, hm = viz(data_mesh;showfacets=false, color=data_viz,colorscheme=:jet1)
wgl.Colorbar(fig3[1,2], colorrange = (minimum(data_viz),maximum(data_viz)), colormap = :jet1, label="v0 (m/s)");
wgl.save("v0_1x1b_$(M)DOFs_$(Int(freq))Hz.html", fig3)

v_abs = sqrt.(v[0M+1:1M].^2 + v[1M+1:2M].^2 + v[2M+1:3M].^2);
data_mesh,data_viz = create_vizualization_data(mesh,abs.(v_abs));
fig4, ax, hm = viz(data_mesh;showfacets=false, color=data_viz,colorscheme=:jet1)
wgl.Colorbar(fig4[1,2], colorrange = (minimum(data_viz),maximum(data_viz)), colormap = :jet1, label="|Vv| (m/s)");
wgl.save("Vv_abs_1x1b_$(M)DOFs_$(Int(freq))Hz.html", fig4)

v_abs = sqrt.(v[0M+1:1M].^2 + v[1M+1:2M].^2 + v[2M+1:3M].^2);
log_vv = log10.(abs.(v_abs));
data_mesh,data_viz = create_vizualization_data(mesh,log_vv);
fig5, ax, hm = viz(data_mesh;showfacets=true, color=data_viz,colorscheme=:jet1)
wgl.Colorbar(fig5[1,2], colorrange = (minimum(data_viz),maximum(data_viz)), colormap = :jet1, label="log10(Vv)");
wgl.save("log10Vv_abs_1x1b_$(M)DOFs_$(Int(freq))Hz.html", fig5)

p_h = -LGM.tau_a/LGM.tau_h.*pa;
log_ph = log10.(abs.(p_h));
data_mesh,data_viz = create_vizualization_data(mesh,log_ph);
fig6, ax, hm = viz(data_mesh;showfacets=false, color=data_viz,colorscheme=:jet1)
wgl.Colorbar(fig6[1,2], colorrange = (minimum(data_viz),maximum(data_viz)), colormap = :jet1, label="log10(Ph)");
wgl.save("log10ph_1x1b_$(M)DOFs_$(Int(freq))Hz.html", fig6)

data_mesh,data_viz = create_vizualization_data(mesh,abs.(v_r));
fig7, ax, hm = viz(data_mesh;showfacets=true, color=data_viz,colorbar=true)
wgl.Colorbar(fig7[1,2], colorrange = (minimum(data_viz),maximum(data_viz)), label="|Vv_r| (m/s)");
wgl.save("Vvr_abs_1x1b_$(M)DOFs_$(Int(freq))Hz.html", fig7)

data_mesh,data_viz = create_vizualization_data(mesh,abs.(v_theta));
fig8, ax, hm = viz(data_mesh;showfacets=true, color=data_viz,colorbar=true)
wgl.Colorbar(fig8[1,2], colorrange = (minimum(data_viz),maximum(data_viz)), label="|Vv_theta| (m/s)");
wgl.save("Vvtheta_abs_1x1b_$(M)DOFs_$(Int(freq))Hz.html", fig8) =#