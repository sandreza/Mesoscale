using JLD2, LinearAlgebra, Statistics
using GLMakie

prefix = "fluxernathy_tracers_restarted_smooth_forcing_j10_k10_averages.jld2"
include("utils.jl")
jl_file = jldopen(prefix, "r+")



time_index = 22
b = get_field("b", time_index, jl_file)
y, z = get_grid(b, jl_file; ghost = 3)

y̅ = avg(y)
z̅ = avg(z)
bz_p = ∂z(b, z)
by_p = ∂y(b, y)

bz = avgy(bz_p)
by = avgz(by_p)

v = get_field("v", time_index, jl_file)
w = get_field("w", time_index, jl_file)

fig = Figure()
options = (; xlabel = "Buoyancy [m/s²]",
    xlabelcolor = :black, ylabel = "Depth [m]",
    ylabelcolor = :black, xlabelsize = 40, ylabelsize = 40,
    xticklabelsize = 25, yticklabelsize = 25,
    xtickcolor = :black, ytickcolor = :black,
    xticklabelcolor = :black, yticklabelcolor = :black,
    titlesize = 50)

indexlist = [50, 100, 150, 200]
ax1 = Axis(fig[1, 1]; options...)
ax2 = Axis(fig[1, 2]; options...)
ax3 = Axis(fig[2, 1]; options...)
ax4 = Axis(fig[2, 2]; options...)

lines!(ax1, b[indexlist[1], :], z, color = :black)
lines!(ax2, b[indexlist[2], :], z, color = :black)
lines!(ax3, b[indexlist[3], :], z, color = :black)
lines!(ax4, b[indexlist[4], :], z, color = :black)



fig = Figure()
options = (; xlabel = "Stratification [1/s²]",
    xlabelcolor = :black, ylabel = "Depth [m]",
    ylabelcolor = :black, xlabelsize = 40, ylabelsize = 40,
    xticklabelsize = 25, yticklabelsize = 25,
    xtickcolor = :black, ytickcolor = :black,
    xticklabelcolor = :black, yticklabelcolor = :black,
    titlesize = 50)

indexlist = [50, 100, 150, 200]
ax1 = Axis(fig[1, 1]; options...)
ax2 = Axis(fig[1, 2]; options...)
ax3 = Axis(fig[2, 1]; options...)
ax4 = Axis(fig[2, 2]; options...)

lines!(ax1, bz_p[indexlist[1], :], z̅, color = :black)
lines!(ax2, bz_p[indexlist[2], :], z̅, color = :black)
lines!(ax3, bz_p[indexlist[3], :], z̅, color = :black)
lines!(ax4, bz_p[indexlist[4], :], z̅, color = :black)

#=
scene, layout  = layoutscene()

lscene = layout[1, 1] = Axis(scene, xlabel = "South to North [m]", 
        xlabelcolor = :black, ylabel = "Depth [m]", 
        ylabelcolor = :black, xlabelsize = 40, ylabelsize = 40,
        xticklabelsize = 25, yticklabelsize = 25,
        xtickcolor = :black, ytickcolor = :black,
        xticklabelcolor  = :black, yticklabelcolor = :black,
        titlesize = 50) 
statestring = "ν ≈  f v'b' / ∂ᶻb / ∂ᶻu"
fieldindex = findall(x->x==statestring, statenames)
u = log10.(abs.(states[fieldindex[1]]))
b = states[4]
statenames[4]

xlims = (0, 1e6)
ylims = (-3000, 0)
xlims = Array(range(xlims[1], xlims[2], length = 4)) 
ylims = Array(range(ylims[1], ylims[2], length = 4)) 

state = u .* 1.0
cmap_rgb = to_colormap(:balance);
clims = quantile.(Ref(state[:]), (0.2,0.8))
heatmap1 = heatmap!(lscene, xlims, ylims, state, interpolate = true, colormap = cmap_rgb, colorrange = clims)
contour!(lscene,0..1e6,-3000..0, b, levels = 40, linewidth = 4, color = :black, alpha = 0.5)
display(scene)
=#