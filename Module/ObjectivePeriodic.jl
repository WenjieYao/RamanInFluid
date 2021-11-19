function MatrixGk(x::Vector, pth; O_mat, kx_s, phys, control, gridap)
    y=zeros(eltype(x),length(x))
    B_mat = MatrixB(pth; control, gridap)
    for ki=1:length(kx_s)
        kx = kx_s[ki]
        kb = VectorValue(kx, 0.0)    
        physk = PhysicalParameters(phys.k, kb, phys.ω, phys.ϵ1, phys.ϵ2, phys.ϵ3, phys.ϵd, phys.μ, phys.R, phys.σs, phys.dpml, phys.LHp, phys.LHn, phys.wg_center, phys.wg_size)
        A_matk = MatrixA(pth; phys=physk, control, gridap)
        y += A_matk' \ (O_mat * (A_matk \ (B_mat * x)))
    end
    return y/length(kx_s)
end

# Parallel sum k
function g_pWk(pW::Vector, kx; O_mat, phys, control, gridap)
    N = num_free_dofs(gridap.FE_U)
    @assert length(pW) == (gridap.np + 2 * N * control.K)
    p0 = zeros(gridap.np)
    for i = 1 : gridap.np
        p0[i] = pW[i]
    end
    W_mat = reinterpret(ComplexF64, reshape(pW[gridap.np + 1 : end], (2 * N, control.K)))

    kb = VectorValue(kx, 0.0)    
    physk = PhysicalParameters(phys.k, kb, phys.ω, phys.ϵ1, phys.ϵ2, phys.ϵ3, phys.ϵd, phys.μ, phys.R, phys.σs, phys.dpml, phys.LHp, phys.LHn, phys.wg_center, phys.wg_size)
        
    pf_vec = pf_p0(p0; control, gridap)
    pfh = FEFunction(gridap.FE_Pf, pf_vec)
    pth = (pf -> Threshold(pf; control)) ∘ pfh
    
    A_mat = MatrixA(pth; phys=physk, control, gridap)
    B_mat = MatrixB(pth; control, gridap)

    dgdp, = Zygote.gradient(p -> g_p(p; O_mat, W_mat, phys=physk, control, gridap), p0)
        
    dgdW = reinterpret(Float64, DgdW(A_mat, W_mat, B_mat, O_mat))
    
    g_value = g_p(p0; O_mat, W_mat, phys=physk, control, gridap)

    g_grad = zeros(length(pW) + 1)
    g_grad[1] = g_value
    g_grad[2 : gridap.np + 1] = dgdp
    g_grad[gridap.np + 2 : end] = 2 * dgdW[:]
    
    return g_grad
end


function gpW_sumk(pW::Vector, grad::Vector;ids, kx_s, dkx, O_mat, phys, control, gridap)
    gk = map_parts(ids) do myids
        mykxs = kx_s[myids]
        mygk = map(kx -> g_pWk(pW, kx; O_mat, phys, control, gridap), mykxs)
        return sum(mygk)
    end
    
    gvalue = 0.0
    igrad = 0
    for gki in gk
        igrad += 1
        if igrad == 1
            gvalue = sum(gki) * dkx / 2 / π * control.Amp
        elseif (length(grad)>0)
            grad[igrad-1] = sum(gki) * dkx / 2 / π * control.Amp
        end
    end
    #gvalue = sum(gk)*dkx/2/π*control.Amp
    open("gvalue.txt", "a") do io
        write(io, "$(gvalue/control.Amp) \n")
    end
    # tc = readdlm("tcount.txt", Int64)[1]
    # open("PV/pvalue$(tc).txt", "w") do iop
    #     for i=1:gridap.np
    #         x_temp = pW[i]
    #         write(iop, "$x_temp \n")
    #     end
    # end
    # tc +=1
    # open("tcount.txt", "w") do iop
    #     write(iop, "$tc \n")
    # end
    return gvalue
end

function gpWk_optimize(p_init, L_local, TOL = 1e-4, MAX_ITER = 500; geo_param, phys, control)
    kx_ini = -π / geo_param.L
    dkx = 2 * π / geo_param.L / control.nkx
    kx_end = π / geo_param.L - dkx
    kx_s = range(kx_ini, kx_end;length = control.nkx)
    backend = SequentialBackend()
    parts = get_part_ids(backend, control.nparts)
    prange = PRange(parts, control.nkx)
    ids = map_parts(get_lid_to_gid, prange.partition)

    if (L_local != geo_param.L)
        MeshGenerator(geo_param, "geometry.msh")
    end
    gridap = GridapFE("geometry.msh", 1, 2, ["DirichletEdges", "DirichletNodes"], ["DesignNodes", "DesignEdges"], ["Target"], ["Source"], control.flag_f)

    # Assemble matrices
    N = num_free_dofs(gridap.FE_U)
    O_mat = MatrixOl(phys.k, phys.ϵ1; gridap)
    
    ##################### Optimize #################
    opt = Opt(:LD_MMA, gridap.np + 2 * N * control.K)
    lb = zeros(gridap.np + 2 * N * control.K)
    lb[gridap.np + 1 : end] = - ones(2 * N * control.K) * Inf
    ub = ones(gridap.np + 2 * N * control.K)
    ub[gridap.np + 1 : end] = ones(2 * N * control.K) * Inf
    opt.lower_bounds = lb
    opt.upper_bounds = ub
    opt.ftol_rel = TOL
    opt.maxeval = MAX_ITER
    opt.max_objective = (pW, grad) -> gpW_sumk(pW, grad; ids, kx_s, dkx, O_mat, phys, control, gridap)
    if (length(p_init) == 0)
        pW_initial = readdlm("pW_opt_value.txt", Float64)
        pW_initial = pW_initial[:]
    else
        if (length(p_init) == 1)
            p_initial = ones(gridap.np) * p_init[1]
        else
            p_initial = p_init
        end
        pW_initial = zeros(gridap.np + 2 * N * control.K)
        pW_initial[1 : gridap.np] = p_initial[:]
        pf_vec = pf_p0(p_initial; control, gridap)
        pfh = FEFunction(gridap.FE_Pf, pf_vec)
        pth = (pf -> Threshold(pf; control)) ∘ pfh

        # A_mat = MatrixA(pth; phys, control, gridap)
        # B_mat = MatrixB(pth; control, gridap)
        
        G_ii, W_raw, info = eigsolve(x -> MatrixGk(x, pth; O_mat, kx_s, phys, control, gridap), rand(ComplexF64, N), 2, :LM; krylovdim = 30)
        W_mat = rand(ComplexF64, N, control.K)
        W_mat = rand(ComplexF64, N, control.K)
        for ib = 1 : 2
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