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
α(ξ) = ξ^(1.2)  # Opacity/alpha along the cmap (0 <= ξ <= 1)
cmap_rgb = to_colormap(:thermometer)
A = α.(range(0, 1, length=length(cmap_rgb)))
clims = ( minimum(T[:,:,zind:end]) * 1.0, maximum(T[:,:,zind:end]))
cmap_rgba = RGBAf0.(cmap_rgb, A)

scene, layout = layoutscene(resolution = (1920, 1080),  backgroundcolor=:black)
lscene = layout[1:2, 1] = LScene(scene)   
volumeobj = volume!(lscene, 0..x, 0..y, 0..z, Ti, colorrange=clims, colormap=cmap_rgba,
                     algorithm=:absorption, absorption=10.0f0,
                     showaxis = false,
        )
        #=
cbar = layout[1, 2] = LColorbar(scene, volumeobj, 
                label = "Buoyancy [m/s²]", labelcolor = :white,
                labelsize = 50,
                width = 50, height = 900, absorption=10.0f0, 
                ticklabelcolor = :black, ticklabelsize = 50) 
                =#
lscene2 = layout[1, 2] = LScene(scene)
volumeobj = volume!(lscene2, 0..x, 0..y, 0..z, Ti, colorrange=clims, colormap=cmap_rgba,
                     algorithm=:absorption, absorption=10.0f0,
                     showaxis = false,
        )   
lscene3 = layout[2, 2] = LScene(scene)
volumeobj = volume!(lscene3, 0..x, 0..y, 0..z, Ti, colorrange=clims, colormap=cmap_rgba,
                     algorithm=:absorption, absorption=10.0f0,
                     showaxis = false,
        )            
# LColorbar(scene, colormap=cmap_rgba,  label = "Activity [spikes/sec]")
scene
##
lscene
##
# LText(scene, "A", textsize = 35, font = "Noto Sans Bold", halign = :right)
##
zoom!(scene.children[1], (100, 00, 00), -1.25, false)
θ = -0.005 * 2π
rotate_cam!(scene.children[1], (θ, 0, 0))
θ2 = -0.005 * 2π
rotate_cam!(scene.children[2], (0, θ2, 0))
θ3 = -0.005 * 2π
rotate_cam!(scene.children[3], (0, 0, θ3))

##
record(scene, "oceananigans_makie.gif", 1:length(t), framerate=10) do n
    obs[] = n
    n == 1 && zoom!(scene.children[1], (0, 0, 0), -1.25, false)
    θ = -0.005 * 2π
    rotate_cam!(scene.children[1], (θ, 0, 0))
    θ2 = -0.005 * 2π
    rotate_cam!(scene.children[2], (0, θ2, 0))
    θ3 = -0.005 * 2π
    rotate_cam!(scene.children[3], (0, 0, θ3))
end
#=


time_index = Node(1)
data = @lift get_data($time_index, 3)

function make_movie(scene)
    n_frames = nt
    θ = 2π / 360

    record(scene, "dry_convection_oceananigans_makie.mp4", 1:n_frames, framerate=30) do n
        @info "frame $n/$n_frames"
        time_index[] = n
        n == 1 && zoom!(scene, (0, 0, 0), -1.25, false)
        rotate_cam!(scene, (θ, 0, 0))
    end
end

scene = volume(0..nx, 0..ny, 0..nz, data, colorrange=(0, 400), colormap=cmap_rgba, algorithm=:absorption, absorption=10.0f0,
               backgroundcolor=cmap_rgb[2], show_axis=false, resolution = (1920, 1080))

make_movie(scene)
=#