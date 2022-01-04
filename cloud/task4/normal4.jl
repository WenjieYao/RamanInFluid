using Gridap, Gridap.Geometry, Gridap.Fields
using Gmsh, GridapGmsh
using DelimitedFiles, Interpolations
using LinearAlgebra, SparseArrays, KrylovKit
using ChainRulesCore, Zygote
using PartitionedArrays
using NLopt

import Gridap.CellData: Interpolable
import ChainRulesCore: rrule
import Gmsh: gmsh


main_path = "/home/gridsan/wyao/Research/RamanInFluid/"
include(main_path*"Materials/Materials.jl")
include(main_path*"Module/Mesh_Periodic.jl")
include(main_path*"Module/Helper.jl")
include(main_path*"Module/GridapFE.jl")
include(main_path*"Module/Control.jl")
include(main_path*"Module/Model.jl")
include(main_path*"Module/Objective.jl")

# Geometry parameters of the mesh
L = 150           # Length of the normal region  
hair = 500        # Height of the air region
hs = 300          # Height of the source location in air
ht = 200          # Height of the target location in air
hd = 200          # Height of design domain
hsub = 100        # Height of substrate domain below design domain
dpml = 100        # Thickness of the PML

# Characteristic length (controls the resolution, smaller the finer)
resol = 100        # Number of points per wavelength
l1 = L/resol      # Air
l2 = l1/2.0       # Design domain
l3 = l1           # PML

meshfile = "geometry.msh"
geo_param = PeriodicGeometry(L, hair, hs, ht, hd, hsub, dpml, l1, l2, l3)
MeshGenerator(geo_param, meshfile)

############  Optimization parameters #############
flag_f = true       # Turn on filter
flag_t = true       # Turn on threshold

# Filter and threshold paramters
r = [0.02 * L, 0.02 * L]  # Filter radius
β = 80.0                  # β∈[1,∞], threshold sharpness
η = 0.5                   # η∈[0,1], threshold center

α = 0.0 / (2 * 1000.0)    # Equivalent loss α = 1/2Q

# Number of subspace
K = 20

# Amplify g for NLopt
Amp = 1

# Sum over kx
nkx = 30
nparts = nkx / 2

Bp = true          # Matrix B depend on parameters?
pv = 1

# Foundary constraint parameters
c = 0#resol^4
lw = r[1]
ls = r[1]
ηe = fηe(lw / r[1])
ηd = fηd(lw / r[1])

control = ControllingParameters(flag_f, flag_t, r, β, η, α, nparts, nkx, K, Amp, Bp, pv, c, ηe, ηd)

gridap = GridapFE(meshfile, 1, 2, ["DirichletEdges", "DirichletNodes"], ["DesignNodes", "DesignEdges"], ["Target"], ["Source"], flag_f)

material = "Silver"
n_λ, k_λ = RefractiveIndex(material, main_path)
λ1 = 532
λ2 = 549
nm1 = n_λ(λ1) + 1im * k_λ(λ1)
nm2 = n_λ(λ2) + 1im * k_λ(λ2)
nf = sqrt(1.77)
μ = 1
R = 1e-10
LHp=[Inf, hair + hd]  # Start of PML for x,y > 0
LHn=[Inf, hsub]       # Start of PML for x,y < 0


ω1 = 2 * π / λ1
phys1 = PhysicalParameters(ω1, nf, nm1, nm1, μ, R, dpml, LHp, LHn, hd)
ω2 = 2 * π / λ2
phys2 = PhysicalParameters(ω2, nf, nm2, nm2, μ, R, dpml, LHp, LHn, hd)

p_init_fun(x) = x[2] < 0.4 * hd ? 1 : 0
lc_temp(v) = ∫(v * x->p_init_fun(x))gridap.dΩ
pc_vec = assemble_vector(lc_temp, gridap.FE_P)
p_init = p_extract(pc_vec; gridap)
p_init[p_init .> 0] .= 1

β_list = [5.0, 10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 80.0, 80.0, 80.0]
# β_list = [80.0, 80.0, 80.0, 80.0, 80.0]

g_opt = 0
for bi = 1 : 10
    β = β_list[bi]
    control = ControllingParameters(flag_f, flag_t, r, β, η, α, nparts, nkx, K, Amp, Bp, pv, c, ηe, ηd)

    if bi == 1
        g_opt, p_opt = g0_p_optimize(p_init, 1e-12, 50; phys1, phys2, control, gridap)
    
    else
        g_opt, p_opt = g0_p_optimize([], 1e-12, 50; phys1, phys2, control, gridap)
    end
    if isfile("p_opt.value.txt")
        run(`rm p_opt_value.txt`)
    end
    open("p_opt_value.txt", "w") do iop
        for i = 1 : length(p_opt)
            p_temp = p_opt[i]
            write(iop, "$p_temp \n")
        end
    end
    open("g_opt_value.txt", "a") do io
        write(io, "$g_opt \n")
    end
end