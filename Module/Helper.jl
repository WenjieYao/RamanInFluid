"""Helper File""" 
# Convert piece-wise constant p (design region) to pvec (whole domain)
function p_extend(p0; gridap)
    p_vec = zeros(num_free_dofs(gridap.FE_P))
    pi = 0
    @assert length(gridap.tags) == num_free_dofs(gridap.FE_P)
    for i = 1 : length(gridap.tags)
        if gridap.tags[i] == gridap.design_tag
            pi += 1
            p_vec[i] = p0[pi]
        end
    end
    p_vec
end

# Extract the design region part from a whole vector
function p_extract(p_vec; gridap)
    p0 = zeros(eltype(p_vec), gridap.np)
    pi = 0
    @assert length(p_vec) == length(gridap.tags)
    for i = 1 : length(gridap.tags)
        if gridap.tags[i] == gridap.design_tag
            pi += 1
            p0[pi] = p_vec[i]
        end
    end
    @assert gridap.np == pi
    p0
end

# Gaussian Distribution function with center x0
function GaussianD(x, x0::AbstractArray, δ::AbstractArray)
    n = length(x)
    @assert (n == length(x0)) && (n == length(δ))
    δn = 1.0
    x_δ = 0.0
    for i = 1 : n
        δn *= √(2π) * δ[i]
        x_δ += ((x[i] - x0[i]) / δ[i])^2
    end
    1.0 / δn * exp(- x_δ / 2.0)
end

# Gaussian Distribution function with center x0 in only Y direction
function GaussianY(x, x0::AbstractArray, δ::AbstractArray)
    n = length(x)
    @assert (n == length(x0)) && (n == length(δ))
    δn = √(2π) * δ[2]
    x_δ = ((x[2] - x0[2]) / δ[2])^2
    return abs(x[1] - x0[1]) <= (δ[1]) ? 1.0 / δn * exp(- x_δ / 2.0) : 0.0
end


# Evaluate the number of contributing values (sum/total > cutoff) for a vector
function num_contributing_values(Gvec::Vector, cutoff = 0.99)
    nmv = length(Gvec)
    gsum = sum(abs.(Gvec))
    gtemp = 0
    for i = 1 : length(Gvec)
        gtemp += abs(Gvec[i])
        if (gtemp) > cutoff * gsum
            nmv = i
            break
        end
    end
    return nmv
end


function Interpolated_Initial_Guess(gridap)
    gridap_guess = GridapFE("InitialGuess/geometry.msh", 1, 2, ["DirichletEdges", "DirichletNodes"], ["DesignNodes", "DesignEdges"], ["Target"], [], true)
    pW_temp = readdlm("InitialGuess/pW_opt_value.txt", Float64)
    pW_temp = pW_temp[:]
    p_init_guess = pW_temp[1 : gridap_guess.np]
    ph_guess = FEFunction(gridap_guess.FE_P, p_extend(p_init_guess; gridap=gridap_guess))
    pfh_init = interpolate(Interpolable(ph_guess),gridap.FE_Pf)
    ph_init = interpolate(Interpolable(pfh_init),gridap.FE_P)
    p_init = p_extract(get_free_dof_values(ph_init); gridap)
    return p_init
end


