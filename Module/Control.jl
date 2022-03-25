"""
This file contains filter, threshold, volume constraint and foundary contraint functions. 
"""
struct ControllingParameters
    flag_f::Bool             # Enable filter
    flag_t::Bool             # Enable threshold
    r::NTuple{2, Float64}    # r = (rx, ry) filter radius
    β::Float64               # Threshold steepness
    η::Float64               # Threshold value
    α::Float64               # Equivalent loss term
    nparts::Int64            # Number of parts for paralell computing
    nkx::Int64               # Number of k points for k integral
    K::Int64                 # Number of contributing eigenvalues
    Amp::Float64             # Tuning amplitude for optimization
    Bp::Bool                 # Matrix B depend on design parameter or not
    pv::Float64              # Volume constraint on design parameter p 
    c::Float64               # Foundary constraint parameter
    ηe::Float64              # Foundary constraint parameter
    ηd::Float64              # Foundary constraint parameter
    hrd::NTuple{2, Float64}  # Height of Raman region
end
# pf_vec = Filter(p_vec)
#a_f(r,u,v) = ∇(v)⊙(TensorValue(r[1]^2,0,0,r[2]^2)⋅∇(u))+v⊙u
a_f(r, u, v) = ∇(v) ⊙ (TensorValue(r[1]^2, 0, 0, r[2]^2) ⋅ ∇(u))

function Filter(p_vec; control, gridap)
    if (control.flag_f)
        ph = FEFunction(gridap.FE_P, p_vec)
        op = AffineFEOperator(gridap.FE_Pf, gridap.FE_Qf) do u, v
            ∫(a_f(control.r, u, v))gridap.dΩ_d + ∫(v * u)gridap.dΩ, ∫(v * ph)gridap.dΩ#, ∫( 0*v )gridap.dΓ_d
          end
        pfh = solve(op)
        return get_free_dof_values(pfh)
    else
        return p_vec
    end
end

# Threshold function
#Threshold(pf;control) = control.flag_t==false ? pf : ((tanh(control.β*control.η)+tanh(control.β*(pf-control.η)))/(tanh(control.β*control.η)+tanh(control.β*(1.0-control.η))))
function Threshold(pfh; control)
    if control.flag_t
        return ((tanh(control.β * control.η) + tanh(control.β * (pfh - control.η)))/
                (tanh(control.β * control.η) + tanh(control.β * (1.0 - control.η))))
    else
        return pfh
    end
end


function VolumeConstraint(pW::Vector, grad::Vector; control, gridap)
    p0 = pW[1 : gridap.np]
    pf_vec = pf_p0(p0; control, gridap)
    pfh = FEFunction(gridap.FE_Pf, pf_vec)
    
    if length(grad) > 0
        grad[gridap.np + 1 : end] = zeros(length(pW) - gridap.np)
        l_temp(v) = ∫(v * ((pf -> Dptdpf(pf; control)) ∘ pfh))gridap.dΩ_d
        grad0 = assemble_vector(l_temp, gridap.FE_Pf)
        grad[1 : gridap.np] = Dgdp(grad0; control, gridap)
    end
    
    ph = (pf -> Threshold(pf; control)) ∘ pfh
    sum(∫(ph)gridap.dΩ_d) - sum(∫(control.pv)gridap.dΩ_d)
end

function fηe(x)
    if (x >= 0) && (x <= 1)
        return  0.25 * x^2 + 0.5
    elseif (x > 1) && (x <= 2)
        return -0.25 * x^2 + x
    else
        return 1.0
    end
end

function fηd(x)
    if (x >= 0) && (x <= 1)
        return -0.25 * x^2 + 0.5
    elseif (x > 1) && (x <= 2)
        return  0.25 * x^2 + 1.0 - x
    else
        return 0
    end
end

fg(g, c) = exp(- c * norm(g)^2)
gc_LW(ph, ηe) = ph > ηe ? 0.0 : (ph - ηe)^2
gc_LS(ph, ηd) = ph < ηd ? 0.0 : (ph - ηd)^2

function LWConstraint(pW::Vector, grad::Vector; control, gridap)
    p0 = pW[1 : gridap.np]
    pf_vec = pf_p0(p0; control, gridap)
    pfh = FEFunction(gridap.FE_Pf, pf_vec)
    ph = (pf -> Threshold(pf; control)) ∘ pfh
    if length(grad) > 0
        # grad[gridap.np + 1 : end] = zeros(length(pW) - gridap.np)
        l_temp(v) = ∫(v * ((ph -> gc_LW(ph, control.ηe)) ∘ pfh) 
                    * ((g -> fg(g, control.c)) ∘ ∇(pfh)) 
                    * (((pf -> Dptdpf(pf; control)) ∘ pfh) 
                    + 2 * ph / (pfh - control.ηe)) 
                    - 2 * control.c * ph * ((ph -> gc_LW(ph, control.ηe)) ∘ pfh) 
                    * ((g -> fg(g, control.c)) ∘ ∇(pfh)) * (∇(v) ⋅ ∇(pfh)))gridap.dΩ_d
        grad0 = assemble_vector(l_temp, gridap.FE_Pf)
        grad[1 : gridap.np] = Dgdp(grad0; control, gridap)
    end
    
    
    sum(∫(ph * ((g -> fg(g, control.c)) ∘ ∇(pfh)) * ((ph -> gc_LW(ph, control.ηe)) ∘ pfh))gridap.dΩ_d)
end
