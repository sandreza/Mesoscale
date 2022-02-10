using JLD2
using GLMakie
using Statistics

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

## now for plotting
options = (; xlabel = "South to North [m]", ylabelsize = 32,
    xlabelsize = 32, xgridstyle = :dash, ygridstyle = :dash, xtickalign = 1,
    xticksize = 30, ytickalign = 1, yticksize = 30,
    xticklabelsize = 30, yticklabelsize = 30, titlesize = 40)

axlist = []
fig = Figure(resolution = (3000, 700))
ax01 = fig[1, 1] = Axis(fig; title = "Meridional Gradient [1/s²]", ylabel = "Depth [m]", options...)
clims01 = symmetric_quantiles(by, 0.01, symmetrize = false)

hm01 = heatmap!(ax01, avg(y), z, by, colormap = :thermometer, interpolate = true, colorrange = clims01)
push!(axlist, ax01)
Colorbar(fig[1, 2], hm01, height = Relative(3 / 4), width = 25, ticklabelsize = 30,
    labelsize = 30, ticksize = 25, tickalign = 1,)


ax02 = fig[1, 3] = Axis(fig; title = "Meriodional Flux [m²/s³]", options...)
clims02 = symmetric_quantiles(vpbp, 0.05, symmetrize = false)
hm02 = heatmap!(ax02, y, z, vpbp, colormap = :thermometer, interpolate = true, colorrange = clims02)
push!(axlist, ax02)
Colorbar(fig[1, 4], hm02, height = Relative(3 / 4), width = 25, ticklabelsize = 30,
    labelsize = 30, ticksize = 25, tickalign = 1,)

ax03 = fig[2, 1] = Axis(fig; ylabel = "Depth [m]", title = "Vertical Gradient [1/s²]", options...)
clims03 = symmetric_quantiles(bz, 0.01, symmetrize = false)
hm03 = heatmap!(ax03, y, avg(z), bz, colormap = :thermometer, interpolate = true, colorrange = clims03)
push!(axlist, ax03)
Colorbar(fig[2, 2], hm03, height = Relative(3 / 4), width = 25, ticklabelsize = 30,
    labelsize = 30, ticksize = 25, tickalign = 1)

ax04 = fig[2, 3] = Axis(fig; ylabel = "Depth [m]", title = "Vertical Flux [m²/s³]", options...)
clims04 = symmetric_quantiles(wpbp, 0.01, symmetrize = false)
hm04 = heatmap!(ax04, y, z, wpbp, colormap = :thermometer, interpolate = true, colorrange = clims04)
push!(axlist, ax04)
Colorbar(fig[2, 4], hm04, height = Relative(3 / 4), width = 25, ticklabelsize = 30,
    labelsize = 30, ticksize = 25, tickalign = 1)

blevels = b[end, 1:2:end]
for ax in axlist
    contour!(ax, y, z, b, levels = blevels, color = :black, linewidth = 3)
end
