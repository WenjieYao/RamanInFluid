"""Gridap finite element implementations"""

struct GridapParameters
    FE_V            # Finite element function space V for fields
    FE_U            # Finite element function space U for fields
    FE_Q            # Finite element function space Q for parameter            
    FE_P            # Finite element function space P for paramter
    FE_Qf           # Finite element function space Qf for filtered paramter
    FE_Pf           # Finite element function space Pf for filtered paramter
    np::Int64       # Number of design parameters
    Ω               # Whole mesh domain
    dΩ              # Numerical integration for whole mesh domain
    dΩ_d            # Numerical integration for design domain
    dΩ_r            # Numerical integration for raman domain
    dΓ_d            # Numerical integration for design boundary
    dΓ_t            # Numerical integration for target boudary
    dΓ_s            # Numerical integration for source boundary
    nb              # Normal vectors of the target boundary
    tags            # Tags of all elements
    design_tag      # Tags of elements in design domain
end


function GridapFE(meshfile, order, degree, diritags, neumanntags, targettags, sourcetags, flag_f = true, ramantag=false)
    model = GmshDiscreteModel(meshfile)

    # Test and trial finite element function space
    # Scalar-valued shape functions,
    # but a complex vector of free DOFs to represent the solution.
    # (this automatically leads to complex sparse matrices when assembling)
    reffe = ReferenceFE(lagrangian, Float64, order)
    FE_V = TestFESpace(model, reffe, dirichlet_tags = diritags, vector_type = Vector{ComplexF64})
    FE_U = FE_V

    # Piece-wise constant parameter FE function space
    p_reffe = ReferenceFE(lagrangian, Float64, 0)
    FE_Q = TestFESpace(model, p_reffe, vector_type = Vector{Float64})
    FE_P = FE_Q

    # Filtered parameter FE function space
    if flag_f
        pf_reffe = ReferenceFE(lagrangian, Float64, 1)
        FE_Qf = TestFESpace(model, pf_reffe, vector_type = Vector{Float64})
        FE_Pf = FE_Qf
    else
        FE_Qf = FE_Q
        FE_Pf = FE_P
    end

    ############### Integration domain ################
    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)

    ############### Get the design region ################
    labels = get_face_labeling(model)
    dimension = num_cell_dims(model)
    tags = get_face_tag(labels, dimension)
    design_tag = get_tag_from_name(labels, "Design")
    cellmask_d = get_face_mask(labels, "Design", dimension)

    Ω_d = Triangulation(model, cellmask_d)
    dΩ_d = Measure(Ω_d, degree)
    # Number of cells in design region (number of design parameters)
    np = num_cells(Ω_d)

    # Desgin/Target/Source line tags
    Γ_d = BoundaryTriangulation(model; tags = neumanntags)
    dΓ_d = Measure(Γ_d, degree)


    Γ_t = BoundaryTriangulation(model; tags = targettags)
    dΓ_t = Measure(Γ_t, degree)
    nb = get_normal_vector(Γ_t)

    Γ_s = BoundaryTriangulation(model; tags = sourcetags)
    dΓ_s = Measure(Γ_s, degree)

    # Raman tags
    if ramantag
        raman_tag = get_tag_from_name(labels, "Raman")
        cellmask_r = get_face_mask(labels, "Raman", dimension)

        Ω_r = Triangulation(model, cellmask_r)
        dΩ_r = Measure(Ω_r, degree)
    else
        dΩ_r = dΩ_d
    end
    gridap = GridapParameters(FE_V,FE_U,FE_Q,FE_P,FE_Qf,FE_Pf,np,Ω,dΩ,dΩ_d,dΩ_r,dΓ_d,dΓ_t,dΓ_s,nb,tags,design_tag)
    return gridap
end