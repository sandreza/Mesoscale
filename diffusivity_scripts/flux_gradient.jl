using JLD2
using GLMakie
using Statistics
using LinearAlgebra

include(pwd() * "/oceananigans_scripts/utils.jl")
include(pwd() * "/diffusivity_scripts/utils_file.jl")

cases = []

for i in 1:5
    push!(cases, "trial" * string(i))
end

for i in 1:16
    if i != 8
        push!(cases, "attempt" * string(i))
    end
end

# cs, cys, czs, vcps, wcps, tracer_string_list = flux_gradient_cases(cases)
# b, by, bz, u, v, w, eke, vpbp, wpbp, y, z = physical_fields()

ti = 1 # tracer index
cmax = norm(cs[ti]) ./ sqrt(length(cs[ti]))
∇c = [[cys[ti] ./ cmax] [czs[ti] ./ cmax]]
u⃗c = [[vcps[ti] ./ cmax] [wcps[ti] ./ cmax]]

## now for plotting
qu_value = 0.03
cmap = :balance
options = (; xlabel = "South to North [m]", ylabelsize = 32,
    xlabelsize = 32, xgridstyle = :dash, ygridstyle = :dash, xtickalign = 1,
    xticksize = 30, ytickalign = 1, yticksize = 30,
    xticklabelsize = 30, yticklabelsize = 30, titlesize = 40)

axlist = []
fig = Figure(resolution = (2000, 1000))
ax01 = fig[1, 1] = Axis(fig; title = "Meridional Gradient [1/s²]", ylabel = "Depth [m]", options...)
clims01 = symmetric_quantiles(∇c[1], qu_value, symmetrize = true)

hm01 = heatmap!(ax01, avg(y), z, ∇c[1], colormap = cmap, interpolate = true, colorrange = clims01)
push!(axlist, ax01)
Colorbar(fig[1, 2], hm01, height = Relative(3 / 4), width = 25, ticklabelsize = 30,
    labelsize = 30, ticksize = 25, tickalign = 1,)

ax02 = fig[1, 3] = Axis(fig; title = "Meriodional Flux [m²/s³]", options...)
clims02 = symmetric_quantiles(u⃗c[1], qu_value, symmetrize = true)
hm02 = heatmap!(ax02, y, z, u⃗c[1], colormap = cmap, interpolate = true, colorrange = clims02)
push!(axlist, ax02)
Colorbar(fig[1, 4], hm02, height = Relative(3 / 4), width = 25, ticklabelsize = 30,
    labelsize = 30, ticksize = 25, tickalign = 1,)

ax03 = fig[2, 1] = Axis(fig; ylabel = "Depth [m]", title = "Vertical Gradient [1/s²]", options...)
clims03 = symmetric_quantiles(∇c[2], qu_value, symmetrize = true)
hm03 = heatmap!(ax03, y, avg(z), ∇c[2], colormap = cmap, interpolate = true, colorrange = clims03)
push!(axlist, ax03)
Colorbar(fig[2, 2], hm03, height = Relative(3 / 4), width = 25, ticklabelsize = 30,
    labelsize = 30, ticksize = 25, tickalign = 1)

ax04 = fig[2, 3] = Axis(fig; ylabel = "Depth [m]", title = "Vertical Flux [m²/s³]", options...)
clims04 = symmetric_quantiles(u⃗c[2], qu_value, symmetrize = true)
hm04 = heatmap!(ax04, y, z, u⃗c[2], colormap = cmap, interpolate = true, colorrange = clims04)
push!(axlist, ax04)
Colorbar(fig[2, 4], hm04, height = Relative(3 / 4), width = 25, ticklabelsize = 30,
    labelsize = 30, ticksize = 25, tickalign = 1)

blevels = b[end, 1:2:end]
for ax in axlist
    contour!(ax, y, z, b, levels = blevels, color = :black, linewidth = 1)
end
