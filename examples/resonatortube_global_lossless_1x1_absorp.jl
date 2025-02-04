# Save boundary and fieldpoint solution for pa_ll (1n-direct), pa (1n-iter) and vv(1n-iter)
# Call with julia --project=.  ./examples/resonatortube_global_losses_1x1_absorp.jl --freq 250 --mesh_file "resonator_tube"

using ArgParse
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--freq"
            help = "frequency"
            default = 259.0
            arg_type = Float64
            arg_type = Float64
        "--mesh_file"
            help = "mesh file"
            default = "resonator_tube_finer"
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
 =#

#===========================================================================================
                        BEM matrix assembly and (iterative) solution for the boundary
===========================================================================================#

LGM = LossyGlobalOuter(mesh,freq;fmm_on=false,depth=1,n=3,progress=false,int_ext=1);

# direct lossless
@info "Computing lossless boundary solution"
#(H-GY)p = 1*Gvs
ρ,c,kₚ,kₐ,kₕ,kᵥ,τₐ,τₕ,ϕₐ,ϕₕ,η,μ = visco_thermal_constants(;freq=freq,S=1);
factor = 1im*2*pi*freq*ρ
pa_ll,hist_pa_ll = gmres(LGM.Ha,LGM.Ga*factor*v0_n;verbose=true,log=true);
#pa_ll=LGM.Ha\(LGM.Ga*factor*v0_n);

p0 = 2*10^(-5);
pa_rms_ll = sqrt.(0.5*real.(pa_ll.*conj.(pa_ll)));
SPL_ll = 20*log10.(pa_rms_ll./p0);

#using DelimitedFiles
#writedlm( "coordinates_tube.txt",  mesh.sources, '\t')
#writedlm( "elements_tube.txt",  mesh.physics_topology, '\t')
#writedlm( "SPL_lossless_tube.txt",  SPL_ll, '\t')

# Saving data
jldsave("results1x1b_$(M)DOFs_$(Int(freq))Hz.JLD2", 
pa_ll=pa_ll,
SPL_ll=SPL_ll,
coordinates=mesh.sources,
elements=mesh.physics_topology,
)

#===========================================================================================
                        BEM matrix assembly and (iterative) solution for the field points
===========================================================================================#
@info "Computing field point solution"

using DelimitedFiles
sources = collect(readdlm("field-points_mics.txt",Float64)')

LGM_fp = BoundaryIntegralEquations.LossyGlobalOuter_fp(mesh,sources,freq;depth=1,n=3,progress=true,int_ext=1);

pa_fp_ll = (LGM_fp.Ga*factor*v0_n - LGM_fp.Ha*pa_ll)./abs.(LGM_fp.C0);

pa_rms_fp_ll = sqrt.(0.5*real.(pa_fp_ll.*conj.(pa_fp_ll)));
SPL_fp_ll = 20*log10.(pa_rms_fp_ll./p0);

#===========================================================================================
                        Absorption coefficients
===========================================================================================#
d = 0.01;
alpha_ll = 1-abs((exp(1im*kₐ*d)-pa_fp_ll[2]/pa_fp_ll[1])/(pa_fp_ll[2]/pa_fp_ll[1]-exp(-1im*kₐ*d))*exp(-1im*kₐ))^2;

# Saving data
jldsave("results_fp_1x1b_$(M)DOFs_$(Int(freq))Hz.JLD2", 
pa_fp_ll=pa_fp_ll,
SPL_fp_ll=SPL_fp_ll,
alpha_ll=alpha_ll
)

println("alpha_ll: $(alpha_ll)")