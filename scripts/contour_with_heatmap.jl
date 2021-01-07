scene, layout  = layoutscene()
lscene = layout[2:4, 2:4] = LAxis(scene, xlabel = "South to North [m]", 
        xlabelcolor = :black, ylabel = "Depth [m]", 
        ylabelcolor = :black, xlabelsize = 40, ylabelsize = 40,
        xticklabelsize = 25, yticklabelsize = 25,
        xtickcolor = :black, ytickcolor = :black,
        xticklabelcolor  = :black, yticklabelcolor = :black,
        titlesize = 50) 

u = states[1]
b = states[4]
statenames[4]

xlims = (0,1)
ylims = (0,1)
state = b
cmap_rgb = to_colormap(:balance)
clims = extrema(b)
heatmap1 = heatmap!(lscene, xlims, ylims, state, interpolate = true, colormap = cmap_rgb, colorrange = clims)


scene = heatmap(v)
contour!(scene, b, levels = 20, colormap = :greys)
c = contour(b, markersize = 1.0)