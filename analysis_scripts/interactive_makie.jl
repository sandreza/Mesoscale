using NetCDF, Plots, GLMakie, AbstractPlotting
using Printf, Statistics
using ImageTransformations, Colors
using AbstractPlotting.MakieLayout


filename = pwd() * "/New_Weno_20_b.nc"
ncinfo(filename)
filename = pwd() * "/New_Weno_20_b.nc"
ncinfo(filename)
xC = Array(NetCDF.open(filename, "xC"))
yC = Array(NetCDF.open(filename, "yC"))
zC = Array(NetCDF.open(filename, "zC"))
t = Array(NetCDF.open(filename, "time"))
sim_day = t ./ 86400
b = NetCDF.open(filename, "b")

##

field_selector = ["u", "v", "w", "b"] #only b is available
ϕ = NetCDF.open(filename, field_selector[end])
timeindexlength = size(ϕ)[end]


scene, layout = layoutscene(resolution = (1920, 1080),  backgroundcolor=:white)
sl1 = layout[2, 2] = LSlider(scene, range = 1:timeindexlength, startvalue = 1)
zind = 20
obs = sl1.value # make index a "node" to allow for movie making

title = "Buoyancy [m/s²] at t="  #need to use sl1.value to access nodes
tmp = @lift $title * string($(sl1.value))
sl2 = layout[2, 1] = LText(scene, tmp) 

field = @lift Array(ϕ[:,:,zind:end,$obs])
Titmp = ϕ[:,:,zind:end,1] 
zonal_average_field = @lift Array(mean(ϕ[:,:,zind:end,$obs], dims = 1)[1,:,:])
hack(t) = t > 10 ? "" : "  " # spacing is silly (otherwise title moves)
day = @lift "3D Buoyancy Field at day " * hack((t[$obs] .- t[1]) ./ 864000) * @sprintf("%2.2f ", (t[$obs] .- t[1]) ./ 864000)
α(ξ) = ξ^(1.2)  # Opacity/alpha along the cmap (0 <= ξ <= 1)
cmap_rgb = to_colormap(:thermometer)
A = α.(range(0, 1, length=length(cmap_rgb)))
clims = ( minimum(b[:,:,zind:end,1]) * 1.0, maximum(b[:,:,zind:end,1]))
cmap_rgba = RGBAf0.(cmap_rgb, A)

x, y, z = size(Titmp)
lscene = layout[1, 1:2] = LScene(scene)   
volumeobj = volume!(lscene, 0..x, 0..y, 0..z, field, colorrange=clims, colormap=cmap_rgba,
                     algorithm=:absorption, absorption=10.0f0,
                     camera = cam3d!
        )

display(scene)
##

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
    tmp3 = layout[0, 1:4] = LText(scene, day, color = :white, textsize = 50)
    layout[0, 2:3] = LText(scene, " ", color = :white, textsize = 50)
    tmp.padding = (100, 50, 0, 0) # left right bottom top
    tmp2.padding = (0, 100, 0, 0) # left right bottom top
    tmp3.padding = (00, 800, 0, 0) # left right bottom top
display(scene)