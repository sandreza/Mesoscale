using NetCDF, Plots, GLMakie, AbstractPlotting
using ImageTransformations, Colors
using AbstractPlotting.MakieLayout
##
filename = pwd() * "/New_Weno_20_b.nc"
ncinfo(filename)
xC = Array(NetCDF.open(filename, "xC"))
yC = Array(NetCDF.open(filename, "yC"))
zC = Array(NetCDF.open(filename, "zC"))
t = Array(NetCDF.open(filename, "time"))
sim_day = t ./ 86400
b = NetCDF.open(filename, "b")
##

T = b[:,:,:,1]
# AbstractPlotting.available_gradients() # see available colors
cmap = :linear_blue_95_50_c20_n256
cmap = :thermometer
clims = (minimum(T) * 1, maximum(T))
cmapa = RGBAf0.(to_colormap(cmap), 1);
cmapa1 = vcat(fill(RGBAf0(0,0,0,0), 10), cmapa);
cmapa2 = vcat(cmapa, fill(RGBAf0(0,0,0,0), 1));
x, y, z = size(T)
scene = volume(0..x, 0..y, 0..z, T, colorrange=clims, algorithm=:absorption, absorption=10.0f0, colormap=cmapa2, show_axis=false)

##
zind = 20
obs = Node(1) # make index a "node" to allow for movie making
Ti = @lift Array(b[:,:,zind:end,$obs])
Ti_s = @lift Array(mean(b[:,:,zind:end,$obs], dims = 1)[1,:,:])
α(ξ) = ξ^(1.2)  # Opacity/alpha along the cmap (0 <= ξ <= 1)
cmap_rgb = to_colormap(:thermometer)
A = α.(range(0, 1, length=length(cmap_rgb)))
clims = ( minimum(T[:,:,zind:end]) * 1.0, maximum(T[:,:,zind:end]))
cmap_rgba = RGBAf0.(cmap_rgb, A)


scene, layout = layoutscene(resolution = (1920, 1080),  backgroundcolor=:black)

lscene = layout[1:3, 1] = LScene(scene)   
volumeobj = volume!(lscene, 0..x, 0..y, 0..z, Ti, colorrange=clims, colormap=cmap_rgba,
                     algorithm=:absorption, absorption=10.0f0,
                     showaxis = true, 
        )

supertitle = layout[2, 2:4] = LText(scene,  " Zonal Average ", textsize = 50, color = :white)

zonal_ax = layout[3,2:4] = LAxis(scene, xlabel = "South to North [m]", 
                         xlabelcolor = :white, ylabel = "Depth [m]", 
                         ylabelcolor = :white, xlabelsize = 40, ylabelsize = 40,
                         xticklabelsize = 25, yticklabelsize = 25,
                         xtickcolor = :white, ytickcolor = :white,
                         xticklabelcolor  = :white, yticklabelcolor = :white, backgroundcolor = :black)
zonalobj = AbstractPlotting.heatmap!(zonal_ax, yC, zC[zind:end], Ti_s, colorrange = clims, 
                                    colormap=cgrad(:thermometer, categorical=true),     
                                    interpolate=true)          
cbar = layout[1, 3] = LColorbar(scene, zonalobj, 
                label = " "^30, labelcolor = :white,
                labelsize = 40, 
                width = 5*5, height = 400, absorption=10.0f0, 
                ticklabelcolor = :white, ticklabelsize = 50) 
tmp = layout[1,2] = LText(scene, "Buoyancy [m/s²] ", color = :white, textsize = 40, rotation = π/2)
tmp2 = layout[1,4] = LText(scene, " ", color = :white, textsize = 50, rotation = π/2)
tmp3 = layout[0,1] = LText(scene, "3D Buoyancy Field", color = :white, textsize = 50)
tmp.padding = (100, 50, 0, 0) # left right bottom top
tmp2.padding = (0, 50, 0, 0) # left right bottom top
tmp3.padding = (200, 300, -10, 0) # left right bottom top
display(scene)
# slice = layout[2, 2]
##
##     
# LColorbar(scene, colormap=cmap_rgba,  label = "Activity [spikes/sec]")
scene
##
# LText(scene, "A", textsize = 35, font = "Noto Sans Bold", halign = :right)
##
record(scene, "oceananigans_makie.gif", 1:80, framerate=10) do n
    obs[] = n
    n == 1 && zoom!(scene.children[1], (0, 0, 0), -1.25, false)
    θ = -0.005 * 2π
    rotate_cam!(scene.children[1], (θ, 0, 0))
end
##