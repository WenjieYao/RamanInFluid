function psfactor(pth, e1=1, e2=1, e3=1)
    return pth * (1 - pth) / (e1 + (e2 - e1)*pth) / (e1 + (e3 - e1)*pth)
end

function MatrixBs(pth, uh; control, gridap, usat = Inf, damp = 1, e1 = 1, e2 = 1, e3 = 1)
    if control.Bp
        B_mat = assemble_matrix(gridap.FE_U, gridap.FE_V) do u, v
            # ∫((1 - pth) * (conj(∇(v) ⋅ ∇(uh)) * ((∇(u) ⋅ ∇(uh)))) / (1 + SaturationFactor(pth, uh, usat, damp, e1, e2)))gridap.dΩ_d + 
            # ∫((x->fr(x, control.hrd[2], control.hrd[1])) * (conj(∇(v) ⋅ ∇(uh)) * ((∇(u) ⋅ ∇(uh)))) / (1 + SaturationFactor(0, uh, usat, damp, e1, e2)))gridap.dΩ_r
            ∫(psfactor(pth, e1, e2, e3) * (conj(∇(uh)) ⋅ ∇(uh)) * ((∇(u) ⋅ ∇(v))) / (1 + SaturationFactor(pth, uh, usat, damp, e1, e2)))gridap.dΩ_d
        end
    else
        B_mat = assemble_matrix(gridap.FE_U, gridap.FE_V) do u, v
            ∫((x->GaussianD(x, [0, 100], [1,1]))* (conj(∇(uh)) ⋅ ∇(uh) * (∇(u) ⋅ ∇(v))))gridap.dΩ
        end
    end
    return B_mat
end

function g0s_pf(pf_vec; kb1, kb2, phys1, phys2, control, gridap, usat = Inf, damp = 1)
    θ1 = asin(kb1[1]/phys1.nf/phys1.ω)
    θ2 = asin(kb2[1]/phys2.nf/phys2.ω)
    pfh = FEFunction(gridap.FE_Pf, pf_vec)
    pth = (pf -> Threshold(pf; control)) ∘ pfh
    A1_mat = MatrixA(pth, kb1; phys=phys1, control, gridap)
    b1_vec = assemble_vector(v->(∫(v*sqrt(cos(θ1)))gridap.dΓ_s), gridap.FE_V)
    u1_vec = A1_mat\b1_vec
    u1h = FEFunction(gridap.FE_U, u1_vec)
    
    B_mat = MatrixBs(pth, u1h; control, gridap, usat, damp, e1=abs2(phys1.nf^2), e2=abs2(phys1.nm^2), e3=abs2(phys2.nm^2))
    # B_mat = MatrixB(pth, u1fix; control, gridap)
    A2_mat = MatrixA(pth, kb2; phys=phys2, control, gridap)
    o_vec = assemble_vector(v->(∫(v*sqrt(cos(θ2)))gridap.dΓ_t), gridap.FE_V)
    v2_vec = A2_mat'\o_vec
    g_temp = v2_vec' * B_mat * v2_vec
    # g_temp = v2fix_vec' * B_mat * v2fix_vec
    @assert abs(imag(g_temp) / real(g_temp)) < 1e-6
    real(g_temp)
end

function rrule(::typeof(g0s_pf), pf_vec; kb1, kb2, phys1, phys2, control, gridap, usat = Inf, damp = 1)
    function U_pullback(dgdg)
      NO_FIELDS, dgdg * Dg0sdpf(pf_vec; kb1, kb2, phys1, phys2, control, gridap, usat, damp)
    end
    g0s_pf(pf_vec; kb1, kb2, phys1, phys2, control, gridap, usat, damp), U_pullback
end

function MatrixdBs(pth, u1h, v2h; gridap, usat = Inf, damp=1, e1 = 1, e2 = 1, e3 = 1)
    if usat == Inf
        B_mat = assemble_matrix(gridap.FE_U, gridap.FE_V) do u, v
            ∫(psfactor(pth, e1, e2, e3) * (conj(∇(v2h)) ⋅ ∇(v2h)) * ((∇(u) ⋅ ∇(v))))gridap.dΩ_d
        end
    else
        temp_x = SaturationFactor(pth, u1h, usat, damp, e1, e2)
        B_mat = assemble_matrix(gridap.FE_U, gridap.FE_V) do u, v
            ∫(psfactor(pth, e1, e2, e3) * (conj(∇(v2h)) ⋅ ∇(v2h)) * ((∇(u) ⋅ ∇(v))) / (1 + temp_x))gridap.dΩ_d - 
            ∫(psfactor(pth, e1, e2, e3) * (conj(∇(v2h)) ⋅ ∇(v2h)) * (conj(∇(u1h)) ⋅ ∇(u1h)) * SatTemp(temp_x, pth, damp, e1, e2) * (∇(v) ⋅ ∇(u)))gridap.dΩ_d
        end
    end
    return B_mat
end

function dpsfactor(pth, e1=1, e2=1, e3=1)
    tempdp= (1 - pth) / (e1 + (e2 - e1)*pth) / (e1 + (e3 - e1)*pth) + pth * (-1 / (e1 + (e2 - e1)*pth) / (e1 + (e3 - e1)*pth) - (e2 - e1) * (1 - pth) / (e1 + (e2 - e1)*pth)/ (e1 + (e2 - e1)*pth) / (e1 + (e3 - e1)*pth) - (e3 - e1) * (1 - pth) / (e1 + (e3 - e1)*pth)/ (e1 + (e2 - e1)*pth) / (e1 + (e3 - e1)*pth))
    return tempdp
end

function pBsdp(pth, u1h, v2h, usat, damp, e1, e2, e3)
    if usat == Inf
        return dpsfactor(pth, e1, e2, e3) * (conj(∇(v2h)) ⋅ ∇(v2h)) * (conj(∇(u1h)) ⋅ ∇(u1h))
    else
        temp_x = SaturationFactor(pth, u1h, usat, damp, e1, e2)
        part1 =  dpsfactor(pth, e1, e2, e3) * (conj(∇(v2h)) ⋅ ∇(v2h)) * (conj(∇(u1h)) ⋅ ∇(u1h)) / (1 + temp_x)
        part2 = psfactor(pth, e1, e2, e3) * (conj(∇(v2h)) ⋅ ∇(v2h)) * (conj(∇(u1h)) ⋅ ∇(u1h))  * SatTemp(temp_x, pth, damp, e1, e2) * (e2-e1) * abs(conj(∇(u1h)) ⋅ ∇(u1h)) / (e1 + (e2 - e1) * pth)
        return part1 + part2
    end
end

function Dg0sdpf(pf_vec; kb1, kb2, phys1, phys2, control, gridap, usat = Inf, damp = 1)
    θ1 = asin(kb1[1]/phys1.nf/phys1.ω)
    θ2 = asin(kb2[1]/phys2.nf/phys2.ω)
    pfh = FEFunction(gridap.FE_Pf, pf_vec)
    pth = (pf -> Threshold(pf; control)) ∘ pfh
    A1_mat = MatrixA(pth, kb1; phys=phys1, control, gridap)
    b1_vec = assemble_vector(v->(∫(v*sqrt(cos(θ1)))gridap.dΓ_s), gridap.FE_V)
    u1_vec = A1_mat \ b1_vec
    u1h = FEFunction(gridap.FE_U, u1_vec)
    
    B_mat = MatrixBs(pth, u1h; control, gridap, usat, damp, e1=abs2(phys1.nf^2), e2=abs2(phys1.nm^2), e3=abs2(phys2.nm^2))
    # B_mat = MatrixB(pth, u1fix; control, gridap)
    A2_mat = MatrixA(pth, kb2; phys=phys2, control, gridap)
    o_vec = assemble_vector(v->(∫(v*sqrt(cos(θ2)))gridap.dΓ_t), gridap.FE_V)
    v2_vec = A2_mat' \ o_vec
    v2h = FEFunction(gridap.FE_U, v2_vec)
    # v2h = FEFunction(gridap.FE_U, v2fix_vec)
    
    v2conjh = FEFunction(gridap.FE_U, conj(v2_vec))
    w2_vec =  A2_mat \ (B_mat * v2_vec)
    w2h = FEFunction(gridap.FE_U, w2_vec) 
    
    if control.Bp
        B_temp = MatrixdBs(pth, u1h, v2h; gridap, usat, damp, e1=abs2(phys1.nf^2), e2=abs2(phys1.nm^2), e3=abs2(phys2.nm^2))
    else
        B_temp = MatrixBs(pth, v2h; control, gridap)
    end
    w1_vec = A1_mat' \ (B_temp * u1_vec)
    w1conjh = FEFunction(gridap.FE_U, conj(w1_vec))
    l_temp(dp) = ∫(real(((pf->Dptdpf(pf; control))∘pfh)*pBsdp(pth, u1h, v2h, usat, damp, abs2(phys1.nf^2), abs2(phys1.nm^2), abs2(phys2.nm^2)) * control.Bp
                         - 2 * 1 * DAdpf(u1h, w1conjh, pfh, kb1; phys=phys1, control)
                         - 2 * 1 * DAdpf(w2h, v2conjh, pfh, kb2; phys=phys2, control)) * dp)gridap.dΩ_d
    dg0dpf = assemble_vector(l_temp, gridap.FE_Pf)
    return dg0dpf
end

# Final objective function
function g0s_p(p0::Vector; kb1, kb2, phys1, phys2, control, gridap, usat = Inf, damp = 1)
    pf_vec = pf_p0(p0; control, gridap)
    g0s_pf(pf_vec; kb1, kb2, phys1, phys2, control, gridap, usat, damp)
end

function g0s_p(p0::Vector, grad::Vector; kb1, kb2, phys1, phys2, control, gridap, usat = Inf, damp = 1)
    if length(grad) > 0
        dgdp, = Zygote.gradient(p -> g0s_p(p; kb1, kb2, phys1, phys2, control, gridap, usat, damp), p0)
        grad[:] = dgdp * control.Amp
    end
    g_value = g0s_p(p0; kb1, kb2, phys1, phys2, control, gridap, usat, damp)
    open("gvalue.txt", "a") do io
        write(io, "$g_value \n")
    end
    return g_value * control.Amp
end

function g0s_p_optimize(p_init, TOL = 1e-4, MAX_ITER = 500, kb1=0, kb2 = 0; phys1, phys2, control, gridap, usat = Inf, damp = 1)
    ##################### Optimize #################
    opt = Opt(:LD_MMA, gridap.np)
    lb = zeros(gridap.np)
    ub = ones(gridap.np)
    opt.lower_bounds = lb
    opt.upper_bounds = ub
    opt.ftol_rel = TOL
    opt.maxeval = MAX_ITER
    opt.max_objective = (p0, grad) -> g0s_p(p0, grad; kb1, kb2, phys1, phys2, control, gridap, usat, damp)
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
        inequality_constraint!(opt, (x, g) -> LSConstraint(x, g; control, gridap), 1e-2)
        inequality_constraint!(opt, (x, g) -> LWConstraint(x, g; control, gridap), 1e-2)
    end

    (g_opt, p_opt, ret) = optimize(opt, p_initial)
    @show numevals = opt.numevals # the number of function evaluations
    
    return g_opt / control.Amp, p_opt
end