scene, layout  = layoutscene()

lscene = layout[1, 1] = Axis(scene, xlabel = "South to North [m]", 
        xlabelcolor = :black, ylabel = "Depth [m]", 
        ylabelcolor = :black, xlabelsize = 40, ylabelsize = 40,
        xticklabelsize = 25, yticklabelsize = 25,
        xtickcolor = :black, ytickcolor = :black,
        xticklabelcolor  = :black, yticklabelcolor = :black,
        titlesize = 50) 
u = states[1]
b = states[4]
statenames[4]

xlims = (0, 1e6)
ylims = (-3000, 0)
xlims = Array(range(xlims[1], xlims[2], length = 4)) 
ylims = Array(range(ylims[1], ylims[2], length = 4)) 

state = b .* 1.0
cmap_rgb = to_colormap(:balance);
clims = extrema(b)
heatmap1 = heatmap!(lscene, xlims, ylims, state, interpolate = true, colormap = cmap_rgb, colorrange = clims)
contour!(lscene,0..1e6,-3000..0, b, levels = 20, linewidth = 4, color = :black, alpha = 0.5)
display(scene)
##
ii = 11
state = states[ii]
b = states[4]
scene = Figure()
cmap_rgb = to_colormap(:viridis);
clims = (quantile(state[:], 0.05) , quantile(state[:], 0.95))
scene = heatmap(state, interpolate = true, colormap = cmap_rgb, colorrange = clims)
cp = contour!(scene, b, levels = 30, linewidth = 4, color = :black, alpha = 0.5)
axis = scene[Axis]
axis.names.axisnames = ["South to North [m]", "Depth [m]"]
newscene = title(scene, statenames[ii] * " with " * statenames[4] * " contours")

ls = colorlegend(scene[end-1], raw = true, camera = campixel!)
scene_final = vbox(newscene, ls)
display(scene_final)