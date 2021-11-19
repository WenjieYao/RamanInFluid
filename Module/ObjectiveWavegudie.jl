function VectorQ(Ey_eig, Mode_norm, gridap)
    l_temp(v) = ∫(v * Ey_eig / Mode_norm) * gridap.dΓ_t[1]
    q_vec = assemble_vector(l_temp, gridap.FE_V)
    return q_vec
end

function g_v(v_vec; B_mat)
    return real(v_vec' * B_mat * v_vec)
end

#v_vec = v_pf(pf)
function v_pf(pf_vec; q_vec, phys, control, gridap)
    pfh = FEFunction(gridap.FE_Pf, pf_vec)
    pth = (pf -> Threshold(pf; control)) ∘ pfh
    A_mat = MatrixA(pth; phys, control, gridap)
    v_vec = A_mat' \ q_vec
    v_vec
end

# Chain Rule : dg/dp = dg/dg*dg/dv*dv/dpf*dpf/dp
# dg/dv=dg/dg*dg/dv
function rrule(::typeof(g_v),v_vec; B_mat)
  function g_pullback(dgdg)
    NO_FIELDS, dgdg * (B_mat * v_vec)
  end
  g_v(v_vec; B_mat), g_pullback
end


# dg/dpf=dg/dv*dv/dpf
function rrule(::typeof(v_pf), pf_vec; q_vec, phys, control, gridap)
    v_vec = v_pf(pf_vec; q_vec, phys, control, gridap)
    function U_pullback(dgdv)
      NO_FIELDS, Dgvdpf(dgdv, v_vec, pf_vec; phys, control, gridap)
    end
    v_vec, U_pullback
end
  
dDpWG(pfh,v2h,dp;control) = control.Dp*real(((pf->dptdpf(pf;control))∘pfh)*v2h*dp)


function Dgvdpf(dgdv, v_vec, pf_vec; phys, control, gridap)
    pfh = FEFunction(gridap.FE_Pf, pf_vec)
    pth = (pf -> Threshold(pf; control)) ∘ pfh
    A_mat = MatrixA(pth; phys, control, gridap)
    w_vec = A_mat \ dgdv
    
    vdh = FEFunction(gridap.FE_U, conj(v_vec))
    vh = FEFunction(gridap.FE_U, v_vec)
    wh = FEFunction(gridap.FE_V, w_vec)
    l_temp(dp) = ∫(real(DBdpf(pfh, vdh, vh; control) - 2 * DAdpf(pfh, vdh, wh; phys, control)) * dp)gridap.dΩ_d
    return assemble_vector(l_temp, gridap.FE_Pf)
end

# Final objective function
function gv_p(p0::Vector; q_vec, B_mat, phys, control, gridap)
    pf_vec = pf_p0(p0; control, gridap)
    v_vec = v_pf(pf_vec; q_vec, phys, control, gridap)
    g_v(v_vec; B_mat)
end

function gv_p(p0::Vector, grad::Vector; q_vec, phys, control, gridap)
    pf_vec = pf_p0(p0; control, gridap)
    pfh = FEFunction(gridap.FE_Pf, pf_vec)
    pth = (pf -> Threshold(pf; control)) ∘ pfh

    B_mat = MatrixB(pth; control, gridap)
    if length(grad) > 0
        dgvdp, = Zygote.gradient(p -> gv_p(p; q_vec, B_mat, phys, control, gridap), p0)
        grad[:] = dgvdp * control.Amp
    end
    g_value = gv_p(p0; q_vec, B_mat, phys, control, gridap)

    #@show g_value
    open("gvalue.txt", "a") do io
        write(io, "$g_value \n")
    end
    # tc = readdlm("tcount.txt", Int64)[1]
    # open("PV/pvalue$(tc).txt", "w") do iop
    #     for i=1:gridap.np
    #         x_temp = p0[i]
    #         write(iop, "$x_temp \n")
    #     end
    # end
    # tc +=1
    # open("tcount.txt", "w") do iop
    #     write(iop, "$tc \n")
    # end

    return g_value * control.Amp
end


function gvp_optimize(p_init, q_vec, TOL = 1e-4, MAX_ITER = 500; phys, control, gridap)
    ##################### Optimize #################
    opt = Opt(:LD_MMA, gridap.np)
    lb = zeros(gridap.np)
    ub = ones(gridap.np)
    opt.lower_bounds = lb
    opt.upper_bounds = ub
    opt.ftol_rel = TOL
    opt.maxeval = MAX_ITER
    opt.max_objective = (p0, grad) -> gv_p(p0, grad; q_vec, phys, control, gridap)
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