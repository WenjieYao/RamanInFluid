NO_FIELDS = ZeroTangent()

function MatrixBe(pth; control, gridap)
    if control.Bp
        B_mat = assemble_matrix(gridap.FE_U, gridap.FE_V) do u, v
            ∫((1-pth) * (∇(u) ⋅ ∇(v)))gridap.dΩ_d + 
            ∫((x->fr(x, control.hrd[2], control.hrd[1])) * (∇(u) ⋅ ∇(v)))gridap.dΩ_r
        end
    else
        B_mat = assemble_matrix(gridap.FE_U, gridap.FE_V) do u, v
            ∫((x->fr(x, control.hrd[2], control.hrd[1])) * (∇(u) ⋅ ∇(v)))gridap.dΩ_r
        end
    end
    return B_mat
end

function ge_pf(pf_vec; kb, phys2, control, gridap)
    pfh = FEFunction(gridap.FE_Pf, pf_vec)
    pth = (pf -> Threshold(pf; control)) ∘ pfh
    
    B_mat = MatrixBe(pth; control, gridap)
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

# Chain Rule : 
# dg/dpf=dg/dg * dg/dpf
function rrule(::typeof(ge_pf), pf_vec; kb, phys2, control, gridap)
    function U_pullback(dgdg)
      NO_FIELDS, dgdg * Dgedpf(pf_vec; kb, phys2, control, gridap)
    end
    ge_pf(pf_vec; kb, phys2, control, gridap), U_pullback
end

function Dgedpf(pf_vec; kb, phys2, control, gridap)
    pfh = FEFunction(gridap.FE_Pf, pf_vec)
    pth = (pf -> Threshold(pf; control)) ∘ pfh
    
    B_mat = MatrixBe(pth; control, gridap)
    # B_mat = MatrixB(pth, u1fix; control, gridap)
    A2_mat = MatrixA(pth, kb; phys=phys2, control, gridap)
    o_vec = VectorO(1, 1; gridap)
    v2_vec = A2_mat' \ o_vec
    v2h = FEFunction(gridap.FE_U, v2_vec)
    # v2h = FEFunction(gridap.FE_U, v2fix_vec)
    
    v2conjh = FEFunction(gridap.FE_U, conj(v2_vec))
    w2_vec =  A2_mat \ (B_mat * v2_vec)
    w2h = FEFunction(gridap.FE_U, w2_vec) 
    
    l_temp(dp) = ∫(real(-((pf->Dptdpf(pf; control))∘pfh) * (∇(v2conjh) ⋅ ∇(v2h)) * control.Bp
                         - 2 * 1 * DAdpf(w2h, v2conjh, pfh, kb; phys=phys2, control)) * dp)gridap.dΩ_d
    dgedpf = assemble_vector(l_temp, gridap.FE_Pf)
    return dgedpf
end


# Final objective function
function ge_p(p0::Vector; kb, phys2, control, gridap)
    pf_vec = pf_p0(p0; control, gridap)
    ge_pf(pf_vec; kb, phys2, control, gridap)
end

function ge_p(p0::Vector, grad::Vector; kb, phys2, control, gridap)
    if length(grad) > 0
        dgdp, = Zygote.gradient(p -> ge_p(p; kb, phys2, control, gridap), p0)
        grad[:] = dgdp * control.Amp
    end
    g_value = ge_p(p0::Vector; kb, phys2, control, gridap)
    open("gvalue.txt", "a") do io
        write(io, "$g_value \n")
    end
    return g_value * control.Amp
end

function ge_p_optimize(p_init, TOL = 1e-4, MAX_ITER = 500; phys2, control, gridap)
    ##################### Optimize #################
    opt = Opt(:LD_MMA, gridap.np)
    lb = zeros(gridap.np)
    ub = ones(gridap.np)
    opt.lower_bounds = lb
    opt.upper_bounds = ub
    opt.ftol_rel = TOL
    opt.maxeval = MAX_ITER
    opt.max_objective = (p0, grad) -> ge_p(p0, grad; kb=0, phys2, control, gridap)
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