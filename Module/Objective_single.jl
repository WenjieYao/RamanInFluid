NO_FIELDS = ZeroTangent()
function gs_pf(pf_vec; kb, phys1, phys2, control, gridap)
    pfh = FEFunction(gridap.FE_Pf, pf_vec)
    pth = (pf -> Threshold(pf; control)) ∘ pfh
    A1_mat = MatrixA(pth, kb; phys=phys1, control, gridap)
    b1_vec = assemble_vector(v->(∫(v)gridap.dΓ_s), gridap.FE_V)
    u1_vec = A1_mat\b1_vec
    u1h = FEFunction(gridap.FE_U, u1_vec)
    
    braman(v) = ∫((x->GaussianD(x, [0,300], [1,1]))*(∇(v) ⋅ ∇(u1h)))gridap.dΩ
    b_raman = assemble_vector(braman, gridap.FE_V)
    A2_mat = MatrixA(pth, kb; phys=phys2, control, gridap)
    v2_vec = A2_mat\b_raman

    O_mat = MatrixOc(phys2.ω, phys2. nf^2; gridap)
    g_temp = v2_vec' * O_mat * v2_vec
    @assert abs(imag(g_temp) / real(g_temp)) < 1e-6
    real(g_temp)
end


# Chain Rule : 
# dg/dpf=dg/dg * dg/dpf
function rrule(::typeof(gs_pf), pf_vec; kb, phys1, phys2, control, gridap)
    function U_pullback(dgdg)
      NO_FIELDS, dgdg * Dgsdpf(pf_vec; kb, phys1, phys2, control, gridap)
    end
    gs_pf(pf_vec; kb, phys1, phys2, control, gridap), U_pullback
end

function Dgsdpf(pf_vec; kb, phys1, phys2, control, gridap)
    pfh = FEFunction(gridap.FE_Pf, pf_vec)
    pth = (pf -> Threshold(pf; control)) ∘ pfh
    A1_mat = MatrixA(pth, kb; phys=phys1, control, gridap)
    b1_vec = assemble_vector(v->(∫(v)gridap.dΓ_s), gridap.FE_V)
    u1_vec = A1_mat \ b1_vec
    u1h = FEFunction(gridap.FE_U, u1_vec)
    
    braman(v) = ∫((x->GaussianD(x, [0,300], [1,1]))*(∇(v) ⋅ ∇(u1h)))gridap.dΩ
    b_raman = assemble_vector(braman, gridap.FE_V)
    A2_mat = MatrixA(pth, kb; phys=phys2, control, gridap)
    v2_vec = A2_mat \ b_raman

    O_mat = MatrixOc(phys2.ω, phys2. nf^2; gridap)
    v2h = FEFunction(gridap.FE_U, v2_vec)
    # v2h = FEFunction(gridap.FE_U, v2fix_vec)
    
    v2h = FEFunction(gridap.FE_U, v2_vec)
    w2_vec =  A2_mat' \ (O_mat * v2_vec)
    w2conjh = FEFunction(gridap.FE_U, conj(w2_vec)) 
    w2h = FEFunction(gridap.FE_U, (w2_vec)) 
    
    b1raman(v) = ∫((x->GaussianD(x, [0,300], [1,1]))*(∇(v) ⋅ ∇(w2h)))gridap.dΩ
    b1_raman = assemble_vector(b1raman, gridap.FE_V)
    w1_vec = A1_mat' \ b1_raman
    w1conjh = FEFunction(gridap.FE_U, conj(w1_vec))
    l_temp(dp) = ∫(real( - 2 * 1 * DAdpf(u1h, w1conjh, pfh, kb; phys=phys1, control)
                         - 2 * 1 * DAdpf(v2h, w2conjh, pfh, kb; phys=phys2, control)) * dp)gridap.dΩ_d
    dgsdpf = assemble_vector(l_temp, gridap.FE_Pf)
    return dgsdpf
end


# Final objective function
function gs_p(p0::Vector; kb, phys1, phys2, control, gridap)
    pf_vec = pf_p0(p0; control, gridap)
    gs_pf(pf_vec; kb, phys1, phys2, control, gridap)
end

function gs_p(p0::Vector, grad::Vector; kb, phys1, phys2, control, gridap)
    if length(grad) > 0
        dgdp, = Zygote.gradient(p -> gs_p(p; kb, phys1, phys2, control, gridap), p0)
        grad[:] = dgdp * control.Amp
    end
    g_value = gs_p(p0::Vector; kb, phys1, phys2, control, gridap)
    open("gvalue.txt", "a") do io
        write(io, "$g_value \n")
    end
    return g_value * control.Amp
end

function gs_p_optimize(p_init, TOL = 1e-4, MAX_ITER = 500; phys1, phys2, control, gridap)
    ##################### Optimize #################
    opt = Opt(:LD_MMA, gridap.np)
    lb = zeros(gridap.np)
    ub = ones(gridap.np)
    opt.lower_bounds = lb
    opt.upper_bounds = ub
    opt.ftol_rel = TOL
    opt.maxeval = MAX_ITER
    opt.max_objective = (p0, grad) -> gs_p(p0, grad; kb=0, phys1, phys2, control, gridap)
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
    end

    (g_opt, p_opt, ret) = optimize(opt, p_initial)
    @show numevals = opt.numevals # the number of function evaluations
    
    return g_opt / control.Amp, p_opt
end


function MatrixOf(gridap)
    # Assemble the matrix
    return assemble_matrix(gridap.FE_U, gridap.FE_V) do u, v
        ∫( (x->GaussianD(x, [0,300], [1,1])) * u * v )gridap.dΩ
    end
end


function gf_pf(pf_vec; kb, phys1, control, gridap)
    pfh = FEFunction(gridap.FE_Pf, pf_vec)
    pth = (pf -> Threshold(pf; control)) ∘ pfh
    A1_mat = MatrixA(pth, kb; phys=phys1, control, gridap)
    b1_vec = assemble_vector(v->(∫(v)gridap.dΓ_s), gridap.FE_V)
    u1_vec = A1_mat\b1_vec

    O_mat = MatrixOf(gridap)
    g_temp = u1_vec' * O_mat * u1_vec
    @assert abs(imag(g_temp) / real(g_temp)) < 1e-6
    real(g_temp)
end



# Chain Rule : 
# dg/dpf=dg/dg * dg/dpf
function rrule(::typeof(gf_pf), pf_vec; kb, phys1, control, gridap)
    function U_pullback(dgdg)
      NO_FIELDS, dgdg * Dgfdpf(pf_vec; kb, phys1, control, gridap)
    end
    gf_pf(pf_vec; kb, phys1, control, gridap), U_pullback
end

function Dgfdpf(pf_vec; kb, phys1, control, gridap)
    pfh = FEFunction(gridap.FE_Pf, pf_vec)
    pth = (pf -> Threshold(pf; control)) ∘ pfh
    A1_mat = MatrixA(pth, kb; phys=phys1, control, gridap)
    b1_vec = assemble_vector(v->(∫(v)gridap.dΓ_s), gridap.FE_V)
    u1_vec = A1_mat \ b1_vec

    O_mat = MatrixOf(gridap)
    # v2h = FEFunction(gridap.FE_U, v2fix_vec)
    
    u1h = FEFunction(gridap.FE_U, u1_vec)
    w1_vec =  A1_mat' \ (O_mat * u1_vec)
    w1conjh = FEFunction(gridap.FE_U, conj(w1_vec)) 
    
    l_temp(dp) = ∫(real( - 2 * 1 * DAdpf(u1h, w1conjh, pfh, kb; phys=phys1, control)) * dp)gridap.dΩ_d
    dgfdpf = assemble_vector(l_temp, gridap.FE_Pf)
    return dgfdpf
end


# Final objective function
function gf_p(p0::Vector; kb, phys1, control, gridap)
    pf_vec = pf_p0(p0; control, gridap)
    gf_pf(pf_vec; kb, phys1, control, gridap)
end

function gf_p(p0::Vector, grad::Vector; kb, phys1, control, gridap)
    if length(grad) > 0
        dgdp, = Zygote.gradient(p -> gf_p(p; kb, phys1, control, gridap), p0)
        grad[:] = dgdp * control.Amp
    end
    g_value = gf_p(p0::Vector; kb, phys1, control, gridap)
    open("gvalue.txt", "a") do io
        write(io, "$g_value \n")
    end
    return g_value * control.Amp
end

function gf_p_optimize(p_init, TOL = 1e-4, MAX_ITER = 500; phys1, control, gridap)
    ##################### Optimize #################
    opt = Opt(:LD_MMA, gridap.np)
    lb = zeros(gridap.np)
    ub = ones(gridap.np)
    opt.lower_bounds = lb
    opt.upper_bounds = ub
    opt.ftol_rel = TOL
    opt.maxeval = MAX_ITER
    opt.max_objective = (p0, grad) -> gf_p(p0, grad; kb=0, phys1, control, gridap)
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
    end

    (g_opt, p_opt, ret) = optimize(opt, p_initial)
    @show numevals = opt.numevals # the number of function evaluations
    
    return g_opt / control.Amp, p_opt
end