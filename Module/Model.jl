"""This file defines all models"""
struct PhysicalParameters
    ω::Float64                   # Frequency
    nf::Float64                  # Fluid (water)
    nm::ComplexF64               # Metal
    ns::ComplexF64               # Substrate
    μ::Float64                   # Magnetic permeability
    R::Float64                   # Reflection of PML
    dpml::Float64                # Thickness of PML
    LHp::Vector{Float64}         # Start of PML for x, y > 0
    LHn::Vector{Float64}         # Start of PML for x, y < 0
    hd::Float64                  # Height of the design region
end

# PML coordinate streching functions
function s_PML(x; phys)
    σ1 = -3 / 4 * log(phys.R) / phys.dpml
    σ2 = -3 / 4 * log(phys.R) / phys.dpml / min(real(phys.ns), 1e-2)
    σ = x[2]>0 ? σ1 : σ2
    xf = [x[1], x[2]]
    u = @. ifelse(xf > 0 , xf - phys.LHp, - xf - phys.LHn)
    return @. ifelse(u > 0,  1 + (1im * σ / phys.ω) * (u / phys.dpml)^2, $(1.0+0im))
end

function ds_PML(x; phys)
    σ1 = -3 / 4 * log(phys.R) / phys.dpml
    σ2 = -3 / 4 * log(phys.R) / phys.dpml / min(real(phys.ns), 1e-2)
    σ = x[2]>0 ? σ1 : σ2
    xf = [x[1], x[2]]
    u = @. ifelse(xf > 0 , xf - phys.LHp, - xf - phys.LHn)
    ds = @. ifelse(u > 0, (2im * σ / phys.ω) * (1 / phys.dpml)^2 * u, $(0.0+0im))
    return ds.*sign.(xf)
end

struct Λ <: Function
    phys
end

function (Λf::Λ)(x)
    s_x,s_y = s_PML(x; Λf.phys)
    return VectorValue(1/s_x, 1/s_y)
end

Fields.∇(Λf::Λ) = x -> TensorValue{2, 2, ComplexF64}(-(Λf(x)[1])^2 * ds_PML(x; Λf.phys)[1], 0, 0, -(Λf(x)[2])^2 * ds_PML(x; Λf.phys)[2])

function ξ0(x; phys)
    if x[2] < 0
        return 1 / phys.ns^2
    elseif x[2] >=0 && x[2] < phys.hd
        return 1 / phys.nf^2 + 0im
    else 
        return 1.0 + 0im
    end
end

ξd(p, nf, nm, α)= 1 / ((nf + (nm - nf) * p)^2 * (1 + α * 1im)) - 1 / nf^2 # in the design region

a_base(u, v, kb; phys) = (x -> ξ0(x; phys)) * ((∇ .* (Λ(phys) * v) - 1im * kb * v) ⊙ ((Λ(phys) .* ∇(u) + 1im * kb * u))) - (phys.ω^2 * phys.μ * (v * u))

a_design(u, v, pth, kb; phys, control) = ((p -> ξd(p, phys.nf, phys.nm, control.α)) ∘ pth) * ((∇(v) - 1im * kb * v) ⊙ (∇(u) + 1im * kb * u)) - phys.ω^2 * 1im * control.α * phys.μ * (v * u)

a_gram(u, v; phys) = (phys.ω^2 * phys.μ * (v * u))

function MatrixA(pth, kb; phys, control, gridap)
    A_mat = assemble_matrix(gridap.FE_U, gridap.FE_V) do u, v
        ∫(a_base(u, v, kb; phys))gridap.dΩ + ∫(a_design(u, v, pth, kb; phys, control))gridap.dΩ_d
    end
    return lu(A_mat)
end

function MatrixB(pth, uh; control, gridap)
    if control.Bp
        B_mat = assemble_matrix(gridap.FE_U, gridap.FE_V) do u, v
            ∫((1-pth) * ((∇(v) ⋅ ∇(uh)) * conj((∇(u) ⋅ ∇(uh)))))gridap.dΩ_d
        end
    else
        B_mat = assemble_matrix(gridap.FE_U, gridap.FE_V) do u, v
            ∫(u*v)gridap.dΓ_s
        end
    end
    return B_mat
end

function MatrixA0(phys, control, gridap)
    return assemble_matrix(gridap.FE_U, gridap.FE_V) do u, v
        ∫(a_gram(u, v; phys))gridap.dΩ + ∫(phys.ω^2 * 1im * control.α * phys.μ * (v * u))gridap.dΩ_d
    end
end


Sp(u, v, ω, ϵ) = 1im / (4 * ω * ϵ) * (u * ∇(v) - v * ∇(u))

# Objective matrix for emitted power from horizontal line
function MatrixOl(ω, ϵ; gridap)
    # Assemble the matrix
    return assemble_matrix(gridap.FE_U, gridap.FE_V) do u, v
        ∫( Sp(u, v, ω, ϵ) ⋅ VectorValue(0, 1) )gridap.dΓ_t
    end
end

# Objective vector for a waveguide mode overlap
function VectorO(Ey_eig, Mode_norm; gridap)
    l_temp(v) = ∫(v * Ey_eig / Mode_norm)gridap.dΓ_t
    o_vec = assemble_vector(l_temp, gridap.FE_V)
    return o_vec
end
