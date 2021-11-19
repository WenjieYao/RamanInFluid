NO_FIELDS = ZeroTangent()
# Objective trace 
function g_pf(pf_vec; O_mat, W_mat, phys, control, gridap)
    pfh = FEFunction(gridap.FE_Pf, pf_vec)
    pth = (pf -> Threshold(pf; control)) ∘ pfh
    A_mat = MatrixA(pth; phys, control, gridap)
    B_mat = MatrixB(pth; control, gridap)
    U_mat = A_mat \ (B_mat * W_mat)
    g_temp = tr((U_mat' * O_mat * U_mat) / (W_mat' * B_mat * W_mat))
    @assert abs(imag(g_temp) / real(g_temp)) < 1e-6
    real(g_temp)
end

#pf = pf_p0(p0)
function pf_p0(p0; control, gridap)
    p_vec = p_extend(p0; gridap)
    pf_vec = Filter(p_vec; control, gridap)
    pf_vec[pf_vec .< 0] .= 0
    pf_vec[pf_vec .> 1] .= 1.0
    pf_vec
end

function MatrixG(x::Vector; A_mat, B_mat, O_mat)
    A_mat' \ (O_mat * (A_mat \ (B_mat * x)))
end

# Chain Rule : 
# dg/dpf=dg/dg * dg/dpf
function rrule(::typeof(g_pf), pf_vec; O_mat, W_mat, phys, control, gridap)
    function U_pullback(dgdg)
      NO_FIELDS, dgdg * Dgdpf(pf_vec, O_mat, W_mat; phys, control, gridap)
    end
    g_pf(pf_vec; O_mat, W_mat, phys, control, gridap), U_pullback
end

Dptdpf(pf; control) = control.flag_t ? control.β * (1.0 - tanh(control.β * (pf - control.η))^2) / (tanh(control.β * control.η) + tanh(control.β * (1.0 - control.η))) : 1.0

Dξdpf(pf, ϵmin, ϵmax; control)= (ϵmin - ϵmax) / (ϵmin + (ϵmax - ϵmin) * Threshold(pf; control))^2 / (1 + control.α * 1im) * Dptdpf(pf; control)

DAdpf(pfh, u, v; phys, control) = ((pf -> Dξdpf(pf, phys.ϵ1, phys.ϵd; control)) ∘ pfh) * ((∇(v) - 1im *phys.kb * v) ⊙ (∇(u) + 1im *phys.kb * u))

DBdpf(pfh, u, v; control) = ((pf -> Dptdpf(pf; control)) ∘ pfh) * (∇(v) ⊙ ∇(u))

function Dgdpf(pf_vec, O_mat, W_mat; phys, control, gridap)
    pfh = FEFunction(gridap.FE_Pf, pf_vec)
    pth = (pf -> Threshold(pf; control)) ∘ pfh
    A_mat = MatrixA(pth; phys, control, gridap)
    B_mat = MatrixB(pth; control, gridap)
    
    U_mat = A_mat \ (B_mat * W_mat)
    WBW = W_mat' * B_mat * W_mat
    dgdU = (O_mat * U_mat) / WBW
    Z_mat = A_mat' \ dgdU
    #Z_mat = A_mat' \ ((O_mat * U_mat) / WBW)
    Wr_mat = W_mat / WBW
    Wl_mat = W_mat * (dgdU' * U_mat)
    #Wl_mat = W_mat / WBW * (U_mat' * O_mat * U_mat)
    
    dgdpf = zeros(num_free_dofs(gridap.FE_Pf))
    for k_i = 1 : control.K
        uh = FEFunction(gridap.FE_U, U_mat[:, k_i])
        zh = FEFunction(gridap.FE_V, conj(Z_mat[:, k_i]))
        if control.Bp
            wh = FEFunction(gridap.FE_U, W_mat[:, k_i])
            wrh = FEFunction(gridap.FE_U, Wr_mat[:, k_i])
            wlh = FEFunction(gridap.FE_V, conj(Wl_mat[:, k_i]))
            l_temp2(dp) = ∫(real(2 * DBdpf(pfh, wh, zh; control) - 2 * DAdpf(pfh, uh, zh; phys, control) - DBdpf(pfh, wrh, wlh; control)) * dp)gridap.dΩ_d
            dgdpf += assemble_vector(l_temp2, gridap.FE_Pf)
        else
            l_temp1(dp) = ∫(- 2 * real(DAdpf(pfh, uh, zh; phys, control)) * dp)gridap.dΩ_d
            dgdpf += assemble_vector(l_temp1, gridap.FE_Pf)
        end
    end
    return dgdpf
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
function g_p(p0::Vector; O_mat, W_mat, phys, control, gridap)
    pf_vec = pf_p0(p0; control, gridap)
    g_pf(pf_vec; O_mat, W_mat, phys, control, gridap)
end

function g_p(p0::Vector, grad::Vector; O_mat, W_mat, phys, control, gridap)
    if length(grad) > 0
        dgdp, = Zygote.gradient(p -> g_p(p; O_mat, W_mat, phys, control, gridap), p0)
        grad[:] = dgdp
    end
    g_value = g_p(p0::Vector; O_mat, W_mat, phys, control, gridap)
    return g_value
end

# Adding W dependence
function DgdW(A, W, B, O)
    WBW = W' * B * W
    U = A \ (B * W)
    B' * (A' \ (O * (U / WBW))) - (B * W / WBW) * (U' * (O * (U / WBW)))
end

function g_pW(pW::Vector, grad::Vector; O_mat, phys, control, gridap)
    N = num_free_dofs(gridap.FE_U)
    @assert length(pW) == (gridap.np + 2 * N * control.K)
    p0 = zeros(gridap.np)
    for i = 1 : gridap.np
        p0[i] = pW[i]
    end
    W_mat = reinterpret(ComplexF64, reshape(pW[gridap.np + 1 : end], (2 * N, control.K)))

    if length(grad) > 0
        pf_vec = pf_p0(p0; control, gridap)
        pfh = FEFunction(gridap.FE_Pf, pf_vec)
        pth = (pf -> Threshold(pf; control)) ∘ pfh
        
        A_mat = MatrixA(pth; phys, control, gridap)
        B_mat = MatrixB(pth; control, gridap)
        
        dgdp, = Zygote.gradient(p -> g_p(p; O_mat, W_mat, phys, control, gridap), p0)
        grad[1 : gridap.np] = dgdp * control.Amp

        dgdW = reinterpret(Float64, DgdW(A_mat, W_mat, B_mat, O_mat))
        grad[gridap.np + 1 : end] = 2 * control.Amp * dgdW[:]
    end
    g_value = g_p(p0; O_mat, W_mat, phys, control, gridap)
    open("gvalue.txt", "a") do io
        write(io, "$g_value \n")
    end
    return g_value * control.Amp
end


function gpW_optimize(p_init, TOL = 1e-4, MAX_ITER = 500, OptAlg = :LD_MMA, IsQl = false; phys, control, gridap)
    # Assemble matrices
    N = num_free_dofs(gridap.FE_U)
    if IsQl
        O_mat = MatrixOl(phys.k, phys.ϵ1; gridap)
    else
        O_mat = MatrixOc(phys.k, phys.ϵ1; gridap)
    end
    
    ##################### Optimize #################
    opt = Opt(OptAlg, gridap.np + 2 * N * control.K)
    lb = zeros(gridap.np + 2 * N * control.K)
    lb[gridap.np + 1 : end] = - ones(2 * N * control.K) * Inf
    ub = ones(gridap.np + 2 * N * control.K)
    ub[gridap.np + 1 : end] = ones(2 * N * control.K) * Inf
    opt.lower_bounds = lb
    opt.upper_bounds = ub
    opt.ftol_rel = TOL
    opt.maxeval = MAX_ITER
    opt.max_objective = (pW, grad) -> g_pW(pW, grad; O_mat, phys, control, gridap)
    if (length(p_init) == 0)
        pW_initial = readdlm("pW_opt_value.txt", Float64)
        pW_initial = pW_initial[:]
    else
        p_initial = p_init
        pW_initial = zeros(gridap.np + 2 * N * control.K)
        pW_initial[1 : gridap.np] = p_initial[:]
        pf_vec = pf_p0(p_initial; control, gridap)
        pfh = FEFunction(gridap.FE_Pf, pf_vec)
        pth = (pf -> Threshold(pf; control)) ∘ pfh
        
        A_mat = MatrixA(pth; phys, control, gridap)
        B_mat = MatrixB(pth; control, gridap)
        
        G_ii, W_raw, info = eigsolve(x -> MatrixG(x; A_mat, B_mat, O_mat), rand(ComplexF64, N), min(5, control.K), :LM; krylovdim = 30)
        W_mat = rand(ComplexF64, N, control.K)
        for ib = 1 : min(5, control.K)
            W_mat[:, ib] = W_raw[ib]
        end
        W_mat = reinterpret(Float64, W_mat)
        pW_initial[gridap.np + 1 : end] = W_mat[:]
        @show abs(sum(G_ii))
    end
    if control.pv < 1
        inequality_constraint!(opt, (x, g) -> VolumeConstraint(x, g; control, gridap), 1e-2)
    end
    if control.c > 0
        equality_constraint!(opt, (x, g) -> LWConstraint(x, g; control, gridap), 1e-8)
    end
    (g_opt, pW_opt, ret) = optimize(opt, pW_initial)
    @show numevals = opt.numevals # the number of function evaluations
    
    return g_opt / control.Amp, pW_opt
end