"""
using GLMakie, DelimitedFiles, Interpolations

include("Materials/Materials.jl")
material = "Gold"
n_λ, k_λ = RefractiveIndex(material)

N = 200
λ_plot = range(200, stop=1200,length=N)
n_plot = zeros(Complex, N)

for i=1:N
    λi = λ_plot[i]
    n_plot[i] = n_λ[λi] + k_λ[λi]*1im
end

scene, layout = layoutscene(resolution=(1400/2, 900/2))
ax = layout[1,1] = Axis(scene)
lin1 = lines!(ax, λ_plot, real(n_plot))
lin2 = lines!(ax, λ_plot, imag(n_plot))
#ax.yscale = log10
ax.xlabel="λ(nm)"
ax.ylabel="n+ik"
ax.title="Dispersion"
#limits!(ax, 0, 50, 1e-16, 1)
axislegend(ax, [lin1, lin2], ["Real", "Imag"], position = :rt,
    orientation = :vertical)
scene

"""
function RefractiveIndex(material)
    filename = "Materials/"*material*".txt"
    RawData = (open(readdlm, filename))

    λ_raw = RawData[2:end,1] * 1e3
    n_raw = RawData[2:end,2]
    k_raw = RawData[2:end,3]

    knots = (λ_raw,)
    itpr = Interpolations.interpolate(knots, n_raw, Gridded(Linear()))
    itpi = Interpolations.interpolate(knots, k_raw, Gridded(Linear()))
    return itpr, itpi
end
