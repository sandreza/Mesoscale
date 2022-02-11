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
numcases = 7 # square root of number of cases
tmp = reshape(collect(1:numcases^2), (numcases,numcases))
axlist = []
fig = Figure(resolution = (2500, 1500))
for ti in 1:numcases^2
    # ti = 1 # tracer index
    ci = argmin(abs.(tmp .- ti))
    ii = ci[1] - 1
    jj = ci[2] - 1
    println("currently on i=", ii, " and j=", jj)

    cmax = norm(cs[ti]) ./ sqrt(length(cs[ti]))
    ∇c = [[cys[ti] ./ cmax] [czs[ti] ./ cmax]]
    u⃗c = [[vcps[ti] ./ cmax] [wcps[ti] ./ cmax]]

    ## now for plotting
    qu_value = 0.03
    cmap = :balance
    options = (; ylabelsize = 0,
        xlabelsize = 0, xgridstyle = :dash, ygridstyle = :dash, xtickalign = 1,
        xticksize = 0, ytickalign = 1, yticksize = 0,
        xticklabelsize = 0, yticklabelsize = 0)

    ax01 = fig[1+2*ii, 1+2*jj] = Axis(fig; options...)
    clims01 = symmetric_quantiles(∇c[1], qu_value, symmetrize = true)

    hm01 = heatmap!(ax01, avg(y), z, ∇c[1], colormap = cmap, interpolate = true, colorrange = clims01)

    ax02 = fig[1+2*ii, 2+2*jj] = Axis(fig; options...)
    clims02 = symmetric_quantiles(u⃗c[1], qu_value, symmetrize = true)
    hm02 = heatmap!(ax02, y, z, u⃗c[1], colormap = cmap, interpolate = true, colorrange = clims02)

    ax03 = fig[2+2*ii, 1+2*jj] = Axis(fig; options...)
    clims03 = symmetric_quantiles(∇c[2], qu_value, symmetrize = true)
    hm03 = heatmap!(ax03, y, avg(z), ∇c[2], colormap = cmap, interpolate = true, colorrange = clims03)

    ax04 = fig[2+2*ii, 2+2*jj] = Axis(fig; options...)
    clims04 = symmetric_quantiles(u⃗c[2], qu_value, symmetrize = true)
    hm04 = heatmap!(ax04, y, z, u⃗c[2], colormap = cmap, interpolate = true, colorrange = clims04)
end

#=
blevels = b[end, 1:2:end]
for ax in axlist
    contour!(ax, y, z, b, levels = blevels, color = :black, linewidth = 1)
end
=#
