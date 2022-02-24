using Gridap, Gridap.Geometry, Gridap.Fields
using Gmsh, GridapGmsh
using DelimitedFiles, Interpolations
using LinearAlgebra, SparseArrays, KrylovKit
using ChainRulesCore, Zygote
using PartitionedArrays
using NLopt
using FileIO

import Images: Gray
import Gridap.CellData: Interpolable
import ChainRulesCore: rrule
import Gmsh: gmsh


main_path = "/home/gridsan/wyao/Research/RamanInFluid/"
include(main_path*"Materials/Materials.jl")
include(main_path*"Module/Mesh_RecCir.jl")
include(main_path*"Module/Helper.jl")
include(main_path*"Module/GridapFE.jl")
include(main_path*"Module/Control.jl")
include(main_path*"Module/Model.jl")
include(main_path*"Module/Objective.jl")
include(main_path*"Module/Objective_single.jl")

L = 600
h1 = 600
h2 = 200
rd = 100
rs = 10
rt = 150
dpml = 300

res = 150
l1 = L/res
l2 = l1/2

hrd = [0, h1/2]
meshfile = "geometry.msh"
geo_param = RecCirGeometry(L, h1, h2, rt, rd, rs, dpml, l1, l2)
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

Bp = false          # Matrix B depend on parameters?
pv = 1

# Foundary constraint parameters
c = 0#resol^4
lw = r[1]
ls = r[1]
ηe = fηe(lw / r[1])
ηd = fηd(lw / r[1])


control = ControllingParameters(flag_f, flag_t, r, β, η, α, nparts, nkx, K, Amp, Bp, pv, c, ηe, ηd, hrd)

gridap = GridapFE(meshfile, 1, 2, ["DirichletEdges", "DirichletNodes"], ["DesignNodes", "DesignEdges"], ["Target"], ["Source"], flag_f)


material = "Ag"
n_λ, k_λ = RefractiveIndex(material,main_path,true)
λ1 = 532
λ2 = 548
nm1 = n_λ(λ1) + 1im * k_λ(λ1)
nm2 = n_λ(λ2) + 1im * k_λ(λ2)
nf = 1
μ = 1
R = 1e-10
LHp=[L/2, h1+h2]   # Start of PML for x,y > 0
LHn=[L/2, 0.1]       # Start of PML for x,y < 0


ω1 = 2 * π / λ1
phys1 = PhysicalParameters(ω1, nf, nm1, nf, μ, R, dpml, LHp, LHn, 0)
ω2 = 2 * π / λ2
phys2 = PhysicalParameters(ω2, nf, nm2, nf, μ, R, dpml, LHp, LHn, 0)


# specify the path to your local image file
img_path = main_path*"Initial/Rasmus_Raman.png"
img = load(img_path)
data = 1.0 .-Float64.(Gray.(img))
function image_to_function(x, data, Lx, Ly)
    Nx, Ny = size(data)
    xi = Int(round((x[1]/Lx + 0.5) * Nx))
    yi = Ny-Int(round(((x[2]-h1/2)/Ly + 0.5) * Ny))
    if xi > 0 && xi <= Nx && yi > 0 && yi <= Ny
        return data[xi,yi]
    else
        return 0.0
    end
end

binit(v) = ∫(v * x->image_to_function(x, data, 280, 280))gridap.dΩ
pc_vec = assemble_vector(binit, gridap.FE_P)
p_init = p_extract(pc_vec; gridap)
p_init[p_init .< 0.5] .= 0
p_init[p_init .> 0.5] .= 1

β_list = [5.0, 10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 80.0, 80.0]
# β_list = [80.0, 80.0, 80.0, 80.0, 80.0]

g_opt = 0
for bi = 1 : 9
    β = β_list[bi]
    control = ControllingParameters(flag_f, flag_t, r, β, η, α, nparts, nkx, K, Amp, Bp, pv, c, ηe, ηd, hrd)

    if bi == 1
        g_opt, p_opt = gs_p_optimize(p_init, 1e-12, 200; phys1, phys2, control, gridap)
    
    else
        g_opt, p_opt = gs_p_optimize([], 1e-12, 200; phys1, phys2, control, gridap)
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
