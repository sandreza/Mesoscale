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
supertitle = layout[1, 2] = LText(scene,  " Buoyancy [m/s²] ", textsize = 50, color = :white)
cbar = layout[2, 2] = LColorbar(scene, volumeobj, 
                label = "Buoyancy [m/s²]", labelcolor = :white,
                labelsize = 50, 
                width = 5*5, height = 700, absorption=10.0f0, 
                ticklabelcolor = :white, ticklabelsize = 50) 

ax = layout[3,2] = LAxis(scene)
AbstractPlotting.heatmap!(ax, 0..y,0..z, Ti_s, colormap=cgrad(:thermometer, categorical=true), interpolate=true)
display(scene)
# slice = layout[2, 2]
##
##     
# LColorbar(scene, colormap=cmap_rgba,  label = "Activity [spikes/sec]")
scene
##
# LText(scene, "A", textsize = 35, font = "Noto Sans Bold", halign = :right)
##
record(scene, "oceananigans_makie.gif", 1:30, framerate=10) do n
    obs[] = n
    n == 1 && zoom!(scene.children[1], (0, 0, 0), -1.25, false)
    θ = -0.005 * 2π
    rotate_cam!(scene.children[1], (θ, 0, 0))
end
##