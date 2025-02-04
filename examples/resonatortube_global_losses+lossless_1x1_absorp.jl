# Save boundary and fieldpoint solution for pa_ll (1n-direct), pa (1n-iter) and vv(1n-iter)
# Call with julia --project=.  ./examples/resonatortube_global_losses_1x1_absorp.jl --freq 250 --mesh_file "resonator_tube"

using ArgParse
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--freq"
            help = "frequency"
            default = 250.0
            arg_type = Float64
            arg_type = Float64
        "--compute_full_solution"
            help = raw"computing full solution"
            default = true
            arg_type = Bool
        "--mesh_file"
            help = "mesh file"
            default = "resonator_tube"
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

#mesh_files =  ["sphere_1m_coarser","sphere_1m_coarse","sphere_1m","sphere_1m_fine","sphere_1m_finer","sphere_1m_4p5k","sphere_1m_extremely_fine","sphere_1m_35k","sphere_1m_77k"];
mesh_path = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes");

#mesh_file = mesh_files[1]
tri_mesh_file = joinpath(mesh_path,mesh_file);
mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[2],
                                                    physics_order=tri_physics_orders[2],int_ext=1);

                                               
hmax = BoundaryIntegralEquations.get_hmax(mesh);
hmin = BoundaryIntegralEquations.get_hmin(mesh);
ρ,c,kₚ,kₐ,kₕ,kᵥ,τₐ,τₕ,ϕₐ,ϕₕ,η,μ = visco_thermal_constants(;freq=freq,S=1);
lambda = 2*pi/kₚ;
elPerLamb = lambda/hmax;

println("Maximum sized elements per wavelength: $(elPerLamb)")
#==========================================================================================
                            Creating excitation vector
==========================================================================================#
xyzb = mesh.sources;
M  = size(xyzb,2);
u₀ = 0.01; # [m/s]
v0x = zeros(M);
v0x[isapprox.(xyzb[1,:], 0.; atol=1E-10)] .= u₀;
v0 = [v0x; zeros(2M)];

v0_n  =  mesh.normals[1,:].*v0[0M+1:1M] + mesh.normals[2,:].*v0[1M+1:2M] + mesh.normals[3,:].*v0[2M+1:3M]; 
#v0_n = -v0x;

#= # alternative way of applying the bc when using the comsol entities
bc_entities = [0]
bc_element_id = Bool.(sum(bc_entities' .∈ mesh.entities,dims=2))[:]
bc_topology = mesh.physics_topology[:,bc_element_id]
bc_node_id = sort(unique(bc_topology)) # gives the same as findall(x->x==true, isapprox.(xyzb[1,:], 0.; atol=1E-10))

# load results instead of computing them
data_file = "results1x1b_$(M)DOFs_$(Int(freq))Hz.JLD2"
f = jldopen(data_file)
pa=f["pa"];
hist_pa=f["hist_pa"];
v=f["v"];
v_r=f["v_r"];
v_theta=f["v_theta"];
hist_dpa=f["hist_dpa"]; =#

#===========================================================================================
                        BEM matrix assembly and (iterative) solution for the boundary
===========================================================================================#

LGM = LossyGlobalOuter(mesh,freq;fmm_on=false,depth=1,n=3,progress=false,int_ext=1);

@info "Computing RHS"
rhs = LGM.Ga*gmres(LGM.inner,(LGM.Dr*v0 - LGM.Nd'*gmres(LGM.Gv,LGM.Hv*v0;verbose=false));verbose=false);

# iterative
@info "Computing lossy boundary solution"
pa,hist_pa = gmres(LGM,rhs;verbose=false,log=true);

# direct lossless
@info "Computing lossless boundary solution"
#(H-GY)p = 1*Gvs
factor = 1im*2*pi*freq*ρ
#pa_ll,hist_pa_ll = gmres(LGM.Ha,LGM.Ga*factor*v0_n;verbose=true,log=true);
pa_ll=LGM.Ha\(LGM.Ga*factor*v0_n);

p0 = 2*10^(-5);
pa_rms = sqrt.(0.5*real.(pa.*conj.(pa)));
pa_rms_ll = sqrt.(0.5*real.(pa_ll.*conj.(pa_ll)));
SPL = 20*log10.(pa_rms./p0);
SPL_ll = 20*log10.(pa_rms_ll./p0);

#using DelimitedFiles
#writedlm( "coordinates_tube.txt",  mesh.sources, '\t')
#writedlm( "elements_tube.txt",  mesh.physics_topology, '\t')
#writedlm( "SPL_lossless_tube.txt",  SPL_ll, '\t')

ph = -LGM.tau_a/LGM.tau_h.*pa;
log_ph = log10.(abs.(ph));

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
#v_theta = cos.(theta).*cos.(phi).*v_tang[0M+1:1M] + cos.(theta).*sin.(phi).*v_tang[1M+1:2M] - sin.(theta).*v_tang[2M+1:3M];
v_theta = sqrt.(v_tang[0M+1:1M].^2 + v_tang[1M+1:2M].^2 + v_tang[2M+1:3M].^2);

v_abs = sqrt.(v[0M+1:1M].^2 + v[1M+1:2M].^2 + v[2M+1:3M].^2);
log_vv = log10.(abs.(v_abs));

# Saving data
jldsave("results1x1b_$(M)DOFs_$(Int(freq))Hz.JLD2", 
pa=pa,
pa_ll=pa_ll,
SPL=SPL,
SPL_ll=SPL_ll,
ph=ph,
log_ph=log_ph,
hist_pa=hist_pa,
#hist_pa_ll=hist_pa_ll,
v=v,
log_vv=log_vv,
v_r=v_r,
v_theta=v_theta,
hist_dpa=hist_dpa,
coordinates=mesh.sources,
elements=mesh.physics_topology,
)

#===========================================================================================
                        BEM matrix assembly and (iterative) solution for the field points
===========================================================================================#
@info "Computing field point solution"

#dpa_ll = LGM.Ga\(LGM.Ha*pa_ll);
dph = gmres(LGM.Gh,LGM.Hh*ph;verbose=false);
dvv = gmres(LGM.Gv,LGM.Hv*v;verbose=false);

using DelimitedFiles
sources = collect(readdlm("field-points_mics.txt",Float64)')

LGM_fp = BoundaryIntegralEquations.LossyGlobalOuter_fp(mesh,sources,freq;depth=1,n=3,progress=true,int_ext=1);

pa_fp = (LGM_fp.Ga*dpa - LGM_fp.Ha*pa)./abs.(LGM_fp.C0);
pa_fp_ll = (LGM_fp.Ga*factor*v0_n - LGM_fp.Ha*pa_ll)./abs.(LGM_fp.C0);
ph_fp = (LGM_fp.Gh*dph - LGM_fp.Hh*ph)./abs.(LGM_fp.C0);
vv_fp = [(LGM_fp.Gv*dvv[0M+1:1M] - LGM_fp.Hv*v[0M+1:1M])./abs.(LGM_fp.C0);(LGM_fp.Gv*dvv[1M+1:2M] - LGM_fp.Hv*v[1M+1:2M])./abs.(LGM_fp.C0);(LGM_fp.Gv*dvv[2M+1:3M] - LGM_fp.Hv*v[2M+1:3M])./abs.(LGM_fp.C0)];

vv_fp = [(LGM_fp.Gv*dvv[0M+1:1M] - LGM_fp.Hv*v[0M+1:1M]);(LGM_fp.Gv*dvv[1M+1:2M] - LGM_fp.Hv*v[1M+1:2M]);(LGM_fp.Gv*dvv[2M+1:3M] - LGM_fp.Hv*v[2M+1:3M])];
N  = size(sources,2);
vv_abs_fp = sqrt.(vv_fp[0N+1:1N].^2 + vv_fp[1N+1:2N].^2 + vv_fp[2N+1:3N].^2);
log_vv_fp = log10.(abs.(vv_abs_fp));

log_ph_fp = log10.(abs.(ph_fp));

pa_rms_fp = sqrt.(0.5*real.(pa_fp.*conj.(pa_fp)));
SPL_fp = 20*log10.(pa_rms_fp./p0);

pa_rms_fp_ll = sqrt.(0.5*real.(pa_fp_ll.*conj.(pa_fp_ll)));
SPL_fp_ll = 20*log10.(pa_rms_fp_ll./p0);

#===========================================================================================
                        Absorption coefficients
===========================================================================================#
d = 0.01;
#ρ,c,kₚ,kₐ,kₕ,kᵥ,τₐ,τₕ,ϕₐ,ϕₕ,η,μ = visco_thermal_constants(;freq=freq,S=1);
alpha = 1-abs((exp(1im*kₐ*d)-pa_fp[2]/pa_fp[1])/(pa_fp[2]/pa_fp[1]-exp(-1im*kₐ*d))*exp(1im*kₐ))^2;
alpha_ll = 1-abs((exp(1im*kₐ*d)-pa_fp_ll[2]/pa_fp_ll[1])/(pa_fp_ll[2]/pa_fp_ll[1]-exp(-1im*kₐ*d))*exp(-1im*kₐ))^2;

# Saving data
jldsave("results_fp_1x1b_$(M)DOFs_$(Int(freq))Hz.JLD2", 
pa_fp=pa_fp,
pa_fp_ll=pa_fp_ll,
SPL_fp=SPL_fp,
SPL_fp_ll=SPL_fp_ll,
ph_fp=ph_fp,
log_ph_fp=log_ph_fp,
vv_fp=vv_fp,
log_vv_fp=log_vv_fp,
alpha=alpha,
alpha_ll=alpha_ll,
)

println("alpha: $(alpha)")
println("alpha_ll: $(alpha_ll)")

#===========================================================================================
                                Plotting and exporting solutions
===========================================================================================#
#= 
using MeshViz
import WGLMakie as wgl
wgl.set_theme!(resolution=(800, 800));
tri_bc_ents = [0,2] # front of tube, resonator
# Plotting boundary conditions
tri_simple_mesh = create_bc_simple_mesh(mesh,tri_bc_ents,false)
viz(tri_simple_mesh;showfacets=true,alpha=0.5)
tri_simple_bc   = create_bc_simple_mesh(mesh,tri_bc_ents[1])
viz!(tri_simple_bc;showfacets=true,color=:red)
tri_simple_bc   = create_bc_simple_mesh(mesh,tri_bc_ents[2])
viz!(tri_simple_bc;showfacets=true,color=:blue)

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

log_vv = log10.(abs.(v_abs));
data_mesh,data_viz = create_vizualization_data(mesh,log_vv);
fig5, ax, hm = viz(data_mesh;showfacets=true, color=data_viz,colorscheme=:jet1)
wgl.Colorbar(fig5[1,2], colorrange = (minimum(data_viz),maximum(data_viz)), colormap = :jet1, label="log10(Vv)");
wgl.save("log10Vv_abs_1x1b_$(M)DOFs_$(Int(freq))Hz.html", fig5)

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
wgl.save("Vvtheta_abs_1x1b_$(M)DOFs_$(Int(freq))Hz.html", fig8)

using WriteVTK
numElem = size(mesh.physics_topology,2)
cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_TRIANGLE, mesh.physics_topology[:,i]) for i=1:numElem] # array that contains all the cells of the mesh 
vtkfile = vtk_grid("resonatortube_$(M)physDOFs", mesh.sources, cells)
vtkfile["n", VTKPointData()] = mesh.normals
vtkfile["t", VTKPointData()] = mesh.tangents
vtkfile["s", VTKPointData()] = mesh.sangents
vtkfile["SPL", VTKPointData()] = SPL
outfiles = vtk_save(vtkfile) 

using DelimitedFiles
writedlm( "coordinates_metatube.txt",  mesh.sources, '\t')
writedlm( "elements_metatube.txt",  mesh.physics_topology, '\t')
writedlm( "SPL_metatube.txt",  SPL, '\t')
writedlm( "SPL_ll_metatube.txt",  SPL_ll, '\t')
writedlm( "log_ph_metatube.txt",  log_ph, '\t')
writedlm( "log_vv_metatube.txt",  log_vv, '\t')

writedlm( "coordinates_fp.txt",  sources, '\t')
writedlm( "SPL_fp_resonatortube.txt",  SPL_fp, '\t')
writedlm( "SPL_fp_ll_resonatortube.txt",  SPL_fp_ll, '\t')
writedlm( "log_ph_fp_resonatortube.txt",  log_ph_fp, '\t')
writedlm( "log_vv_fp_resonatortube.txt",  log_vv_fp, '\t')
writedlm( "alpha.txt",  alpha, '\t')
writedlm( "alpha_ll.txt",  alpha_ll, '\t') =#