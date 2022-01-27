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