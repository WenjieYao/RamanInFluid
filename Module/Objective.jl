NO_FIELDS = ZeroTangent()
function g0_pf(pf_vec; kb, phys1, phys2, control, gridap)
    pfh = FEFunction(gridap.FE_Pf, pf_vec)
    pth = (pf -> Threshold(pf; control)) ∘ pfh
    A1_mat = MatrixA(pth, kb; phys=phys1, control, gridap)
    b1_vec = assemble_vector(v->(∫(v)gridap.dΓ_s), gridap.FE_V)
    u1_vec = A1_mat\b1_vec
    u1h = FEFunction(gridap.FE_U, u1_vec)
    
    B_mat = MatrixB(pth, u1h; control, gridap)
    # B_mat = MatrixB(pth, u1fix; control, gridap)
    A2_mat = MatrixA(pth, kb; phys=phys2, control, gridap)
    o_vec = VectorO(1, 1; gridap)
    v2_vec = A2_mat'\o_vec
    g_temp = v2_vec' * B_mat * v2_vec
    # g_temp = v2fix_vec' * B_mat * v2fix_vec
    @assert abs(imag(g_temp) / real(g_temp)) < 1e-6
    real(g_temp)
    # v2h = FEFunction(gridap.FE_U, v2_vec)
    # sum(∫((1 - 1 * pth) * abs2(∇(v2h) ⋅ ∇(u1h)))gridap.dΩ_d)
    # sum(∫((1 - pth) * abs2(v2h * u1h))gridap.dΩ_d)
end

#pf = pf_p0(p0)
function pf_p0(p0; control, gridap)
    p_vec = p_extend(p0; gridap)
    pf_vec = Filter(p_vec; control, gridap)
    pf_vec[pf_vec .< 0] .= 0
    pf_vec[pf_vec .> 1] .= 1.0
    pf_vec
end

# Chain Rule : 
# dg/dpf=dg/dg * dg/dpf
function rrule(::typeof(g0_pf), pf_vec; kb, phys1, phys2, control, gridap)
    function U_pullback(dgdg)
      NO_FIELDS, dgdg * Dg0dpf(pf_vec; kb, phys1, phys2, control, gridap)
    end
    g0_pf(pf_vec; kb, phys1, phys2, control, gridap), U_pullback
end

Dptdpf(pf; control) = control.flag_t ? control.β * (1.0 - tanh(control.β * (pf - control.η))^2) / (tanh(control.β * control.η) + tanh(control.β * (1.0 - control.η))) : 1.0

Dξdpf(pf, nf, nm, control)= 2 * (nf - nm) / (nf + (nm - nf) * Threshold(pf; control))^3 / (1 + control.α * 1im) * Dptdpf(pf; control)

DAdpf(u, v, pfh, kb; phys, control) = ((p -> Dξdpf(p, phys.nf, phys.nm, control)) ∘ pfh) * ((∇(v) - 1im * kb * v) ⊙ (∇(u) + 1im * kb * u))

function Dg0dpf(pf_vec; kb, phys1, phys2, control, gridap)
    pfh = FEFunction(gridap.FE_Pf, pf_vec)
    pth = (pf -> Threshold(pf; control)) ∘ pfh
    A1_mat = MatrixA(pth, kb; phys=phys1, control, gridap)
    b1_vec = assemble_vector(v->(∫(v)gridap.dΓ_s), gridap.FE_V)
    u1_vec = A1_mat \ b1_vec
    u1h = FEFunction(gridap.FE_U, u1_vec)
    
    B_mat = MatrixB(pth, u1h; control, gridap)
    # B_mat = MatrixB(pth, u1fix; control, gridap)
    A2_mat = MatrixA(pth, kb; phys=phys2, control, gridap)
    o_vec = VectorO(1, 1; gridap)
    v2_vec = A2_mat' \ o_vec
    v2h = FEFunction(gridap.FE_U, v2_vec)
    # v2h = FEFunction(gridap.FE_U, v2fix_vec)
    
    v2conjh = FEFunction(gridap.FE_U, conj(v2_vec))
    w2_vec =  A2_mat \ (B_mat * v2_vec)
    w2h = FEFunction(gridap.FE_U, w2_vec) 
    
    B_temp = MatrixB(pth, v2h; control, gridap)
    w1_vec = A1_mat' \ (B_temp * u1_vec)
    w1conjh = FEFunction(gridap.FE_U, conj(w1_vec))
    l_temp(dp) = ∫(real(-((pf->Dptdpf(pf; control))∘pfh) * abs2(∇(v2h) ⋅ ∇(u1h)) * control.Bp
                         - 2 * 1 * DAdpf(u1h, w1conjh, pfh, kb; phys=phys1, control)
                         - 2 * 1 * DAdpf(w2h, v2conjh, pfh, kb; phys=phys2, control)) * dp)gridap.dΩ_d
    dg0dpf = assemble_vector(l_temp, gridap.FE_Pf)
    return dg0dpf
end

# dg/dp=dg/dpf*dpf/dp
function rrule(::typeof(pf_p0), p0; control, gridap)
  function pf_pullback(dgdpf)
    NO_FIELDS, Dgdp(dgdpf; control, gridap)
  end
  pf_p0(p0; control, gridap), pf_pullback
end

function Dgdp(dgdpf; control, gridap)
    if control.flag_f
        Af = assemble_matrix(gridap.FE_Pf, gridap.FE_Qf) do u, v
            ∫(a_f(control.r, u, v))gridap.dΩ_d + ∫(v * u)gridap.dΩ
        end
        λvec = Af' \ dgdpf
        λh = FEFunction(gridap.FE_Pf, λvec)
        l_temp(dp) = ∫(λh * dp)gridap.dΩ
        return p_extract(assemble_vector(l_temp, gridap.FE_P); gridap)
    else
        return p_extract(dgdpf; gridap)
    end
end


# Final objective function
function g0_p(p0::Vector; kb, phys1, phys2, control, gridap)
    pf_vec = pf_p0(p0; control, gridap)
    g0_pf(pf_vec; kb, phys1, phys2, control, gridap)
end

function g0_p(p0::Vector, grad::Vector; kb, phys1, phys2, control, gridap)
    if length(grad) > 0
        dgdp, = Zygote.gradient(p -> g0_p(p; kb, phys1, phys2, control, gridap), p0)
        grad[:] = dgdp * control.Amp
    end
    g_value = g0_p(p0::Vector; kb, phys1, phys2, control, gridap)
    open("gvalue.txt", "a") do io
        write(io, "$g_value \n")
    end
    return g_value * control.Amp
end

function g0_p_optimize(p_init, TOL = 1e-4, MAX_ITER = 500; phys1, phys2, control, gridap)
    ##################### Optimize #################
    opt = Opt(:LD_MMA, gridap.np)
    lb = zeros(gridap.np)
    ub = ones(gridap.np)
    opt.lower_bounds = lb
    opt.upper_bounds = ub
    opt.ftol_rel = TOL
    opt.maxeval = MAX_ITER
    opt.max_objective = (p0, grad) -> g0_p(p0, grad; kb=0, phys1, phys2, control, gridap)
    if (length(p_init)==0)
        p_initial = readdlm("p_opt_value.txt", Float64)
        p_initial = p_initial[:]
    else
        p_initial = p_init[:]
    end
    if control.pv < 1
        inequality_constraint!(opt, (x, g) -> VolumeConstraint(x, g; control, gridap), 1e-2)
    end
    if control.c > 0
        equality_constraint!(opt, (x, g) -> LWConstraint(x, g; control, gridap), 1e-8)
    end

    (g_opt, p_opt, ret) = optimize(opt, p_initial)
    @show numevals = opt.numevals # the number of function evaluations
    
    return g_opt / control.Amp, p_opt
end