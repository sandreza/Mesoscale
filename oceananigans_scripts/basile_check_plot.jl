include("basile_check.jl")

using GLMakie

##
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



fig2 = Figure()
options = (; xlabel = "Stratification [1/s²]",
    xlabelcolor = :black, ylabel = "Depth [m]",
    ylabelcolor = :black, xlabelsize = 40, ylabelsize = 40,
    xticklabelsize = 25, yticklabelsize = 25,
    xtickcolor = :black, ytickcolor = :black,
    xticklabelcolor = :black, yticklabelcolor = :black,
    titlesize = 50)

indexlist = [50, 100, 150, 200]
ax1 = Axis(fig2[1, 1]; options...)
ax2 = Axis(fig2[1, 2]; options...)
ax3 = Axis(fig2[2, 1]; options...)
ax4 = Axis(fig2[2, 2]; options...)

lines!(ax1, bz_p[indexlist[1], :], z̅, color = :black)
lines!(ax2, bz_p[indexlist[2], :], z̅, color = :black)
lines!(ax3, bz_p[indexlist[3], :], z̅, color = :black)
lines!(ax4, bz_p[indexlist[4], :], z̅, color = :black)


function mplot(z̅, bz_p; xlabel = "Stratification [1/s²]")
    fig2 = Figure()
    options = (; xlabel = xlabel,
        xlabelcolor = :black, ylabel = "Depth [m]",
        ylabelcolor = :black, xlabelsize = 40, ylabelsize = 40,
        xticklabelsize = 25, yticklabelsize = 25,
        xtickcolor = :black, ytickcolor = :black,
        xticklabelcolor = :black, yticklabelcolor = :black,
        titlesize = 50)

    indexlist = [50, 100, 150, 200]
    ax1 = Axis(fig2[1, 1]; options...)
    ax2 = Axis(fig2[1, 2]; options...)
    ax3 = Axis(fig2[2, 1]; options...)
    ax4 = Axis(fig2[2, 2]; options...)

    lines!(ax1, bz_p[indexlist[1], :], z̅, color = :black)
    lines!(ax2, bz_p[indexlist[2], :], z̅, color = :black)
    lines!(ax3, bz_p[indexlist[3], :], z̅, color = :black)
    lines!(ax4, bz_p[indexlist[4], :], z̅, color = :black)
    display(fig2)
    return fig2
end


##
using Statistics
options = (; xlabel = "barotropic projections",
    xlabelcolor = :black, ylabel = "Depth [m]",
    ylabelcolor = :black, xlabelsize = 40, ylabelsize = 40,
    xticklabelsize = 25, yticklabelsize = 25,
    xtickcolor = :black, ytickcolor = :black,
    xticklabelcolor = :black, yticklabelcolor = :black,
    titlesize = 50)

fig2 = Figure()
ax1 = Axis(fig2[1, 1]; options...)


truth = lines!(ax1, uu, z[zi_b-1:zi_s], color = :black)
b0 = lines!(ax1, û[end] * V[:, end], z[zi_b-1:zi_s], color = :red)
b1 = lines!(ax1, û[end] * V[:, end] + û[end-1] * V[:, end-1], z[zi_b-1:zi_s], color = :blue)
b2 = lines!(ax1, û[end] * V[:, end] + û[end-1] * V[:, end-1] + û[end-2] * V[:, end-2], z[zi_b-1:zi_s], color = :purple)
axislegend(ax1, [truth, b0, b1, b2], ["truth", "b0", "b1", "b2"], position = :rb)

##
options = (; xlabel = "y [m]",
    xlabelcolor = :black, ylabel = "Diffusivity",
    ylabelcolor = :black, xlabelsize = 40, ylabelsize = 40,
    xticklabelsize = 25, yticklabelsize = 25,
    xtickcolor = :black, ytickcolor = :black,
    xticklabelcolor = :black, yticklabelcolor = :black,
    titlesize = 50)

fig2 = Figure()
ax3 = Axis(fig2[1, 1]; options...)

scatter!(ax3, [D̅list[:]...])
ylims!(ax3, (0, 0.25))

##
options = (; xlabel = "y indices",
    xlabelcolor = :black, ylabel = "First Baroclinic Mode Amplitude",
    ylabelcolor = :black, xlabelsize = 40, ylabelsize = 40,
    xticklabelsize = 25, yticklabelsize = 25,
    xtickcolor = :black, ytickcolor = :black,
    xticklabelcolor = :black, yticklabelcolor = :black,
    titlesize = 50)

fig2 = Figure()
ax3 = Axis(fig2[1, 1]; options...)

scatter!(ax3, [𝒰list[:]...])
# ylims!(ax3, (0, 0.25))
##
options = (; xlabel = "y indices",
    xlabelcolor = :black, ylabel = "First Baroclinic Mode Amplitude",
    ylabelcolor = :black, xlabelsize = 40, ylabelsize = 40,
    xticklabelsize = 25, yticklabelsize = 25,
    xtickcolor = :black, ytickcolor = :black,
    xticklabelcolor = :black, yticklabelcolor = :black,
    titlesize = 50)

fig2 = Figure()
ax3 = Axis(fig2[1, 1]; options...)

scatter!(ax3, [λlist[:]...])