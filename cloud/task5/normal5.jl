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

init_ratio = 1.0
init_ratioL = 2.0
init_value = 1.0
init_r = 5
usat = 1e4

material = "Ag"
n_λ, k_λ = RefractiveIndex(material,main_path,true)
λ1 = 532
λ2 = 549
nm1 = n_λ(λ1) + 1im * k_λ(λ1)
nm2 = n_λ(λ2) + 1im * k_λ(λ2)
nf = sqrt(1.77)
μ = 1
R = 1e-10

hr = (λ1+λ2)/nf/4          # Height of Raman molecule
# Geometry parameters of the mesh
L = 400           # Length of the normal region
hair = 500 + hr       # Height of the air region
hs = 300 + hr         # Height of the source location in air
ht = 200 + hr         # Height of the target location in air
hd = 200          # Height of design domain
hsub = 100        # Height of substrate domain below design domain
dpml = 300        # Thickness of the PML
hrd = (hd, hr)
# Characteristic length (controls the resolution, smaller the finer)
resol = 20        # Number of points per wavelength
l1 = 20      # Air
l2 = 1       # Design domain
l3 = l1           # PML

meshfile = "geometry.msh"
geo_param = PeriodicGeometry(L, hair, hs, ht, hd, hsub, dpml, l1, l2, l3)
MeshGenerator(geo_param, meshfile)

LHp=(Inf, hair + hd)  # Start of PML for x,y > 0
LHn=(Inf, hsub)       # Start of PML for x,y < 0


ω1 = 2 * π / λ1
phys1 = PhysicalParameters(ω1, nf, nm1, nm1, μ, R, dpml, LHp, LHn, hd)
ω2 = 2 * π / λ2
phys2 = PhysicalParameters(ω2, nf, nm2, nm2, μ, R, dpml, LHp, LHn, hd)

############  Optimization parameters #############
flag_f = true       # Turn on filter
flag_t = true       # Turn on threshold

# Filter and threshold paramters
r = (init_r, init_r)  # Filter radius
β = 32.0                  # β∈[1,∞], threshold sharpness
η = 0.5                   # η∈[0,1], threshold center

α = 1.0 / (2 * 1000.0)    # Equivalent loss α = 1/2Q

# Number of subspace
K = 20

# Amplify g for NLopt
Amp = 1e-5

# Sum over kx
nkx = 30
nparts = nkx / 2

Bp = true          # Matrix B depend on parameters?
pv = 1

# Foundary constraint parameters
c = (r[1])^4
lw = r[1]
ls = r[1]
ηe = 0.75#fηe(lw / r[1])
ηd = 0.25#fηd(lw / r[1])

control = ControllingParameters(flag_f, flag_t, r, β, η, α, nparts, nkx, K, Amp, Bp, pv, c, ηe, ηd, hrd)

gridap = GridapFE(meshfile, 1, 2, ["DirichletEdges", "DirichletNodes"], ["DesignNodes", "DesignEdges"], ["Target"], ["Source"], flag_f, true)

function p_triangle(x, h, L)
    dx = abs(x[1])
    dy = x[2]
    inner = 0
    if dy <= (1-2*dx/L)*h
        inner = 1
    end
    return inner
end

kb = 0
p_trunc(x, ratio) = x[2] < (ratio * hd) ? 1 : 0
# binitialfunc(v) = ∫(v * x->p_trunc(x, init_ratio))gridap.dΩ
# binitialfunc(v) = ∫(v * x->p_bowtie(x, 20, 80, L, hd))gridap.dΩ
binitialfunc(v) = ∫(v * x->p_triangle(x, init_ratio * hd, init_ratioL * L))gridap.dΩ
pc_vec = assemble_vector(binitialfunc, gridap.FE_P)
p_init = p_extract(pc_vec; gridap)
p_init[p_init .< 0.1] .= 0
p_init[p_init .> 0.1] .= init_value

β_list = [8.0, 8.0, 16.0, 16.0, 32.0, 32.0, 32.0]
Q_list = [10.0, 50.0, 100.0, 500.0, 1000.0, 1000.0, 1000.0]
d_list = [1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 1e-2, 1e-2]

g_opt = 0
for bi = 1 : 6
    β = β_list[bi]
    α = 1/(2*Q_list[bi])
    damp = d_list[bi]
    if bi < 6
        c = 0
        control = ControllingParameters(flag_f, flag_t, r, β, η, α, nparts, nkx, K, Amp, Bp, pv, c, ηe, ηd, hrd)
    else
        c = (r[1])^4
        control = ControllingParameters(flag_f, flag_t, r, β, η, α, nparts, nkx, K, Amp, Bp, pv, c, ηe, ηd, hrd)
    end

    if bi == 1
        g_opt, p_opt = g0_p_optimize(p_init, 1e-12, 100; phys1, phys2, control, gridap, usat, damp)
    
    else
        g_opt, p_opt = g0_p_optimize([], 1e-12, 100; phys1, phys2, control, gridap, usat, damp)
    end
    # if bi == 1
    #     g_opt, p_opt = g0_p_optimize(p_init, 1e-12, 100; phys1, phys2, control, gridap)
    
    # else
    #     g_opt, p_opt = g0_p_optimize([], 1e-12, 100; phys1, phys2, control, gridap)
    # end
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