using NetCDF, Plots, GLMakie, AbstractPlotting
using ImageTransformations, Colors
##
filename = pwd() * "/New_Weno_20_b.nc"
ncinfo(filename)
x = Array(NetCDF.open(filename, "xC"))
y = Array(NetCDF.open(filename, "yC"))
z = Array(NetCDF.open(filename, "zC"))
t = Array(NetCDF.open(filename, "time"))
sim_day = t ./ 86400
z = (z[2:end] + z[1:end-1]) / 2
b = NetCDF.open(filename, "b")
##
T = Array(b[:,:,:,1])
T2 = Array(b[:,:,:,end])
obs = Node(1)
Ti = @lift Array(b[:,:,:,$obs])
# AbstractPlotting.available_gradients() # see available colors
cmap = :linear_blue_95_50_c20_n256
cmap = :thermometer
clims = (minimum(T) * 1, maximum(T))
cmapa = RGBAf0.(to_colormap(cmap), 1);
cmapa1 = vcat(fill(RGBAf0(0,0,0,0), 10), cmapa);
cmapa2 = vcat(cmapa, fill(RGBAf0(0,0,0,0), 1));
x, y, z = size(T)
scene = volume(0..x, 0..y, 0..z, T, colorrange=clims, algorithm=:absorption, absorption=5.0f0, colormap=cmapa2, show_axis=false)

##
α(ξ) = ξ  # Opacity/alpha along the cmap (0 <= ξ <= 1)
cmap_rgb = to_colormap(:thermometer)
A = α.(range(0, 1, length=length(cmap_rgb)))
cmap_rgba = RGBAf0.(cmap_rgb, A)
scene = volume(0..x, 0..y, 0..z, Ti, colorrange=clims, colormap=cmap_rgba, algorithm=:absorption, absorption=10.0f0,
               backgroundcolor=:black, show_axis=false)

##
zoom!(scene, (100, 00, 00), -1.25, false)
θ = 0.015 * 2π
rotate_cam!(scene, (θ, 0, 0))

##
record(scene, "oceananigans_makie.gif", 1:length(t), framerate=30) do n
    obs[] = n
    n == 1 && zoom!(scene, (0, 0, 0), -1.25, false)
    rotate_cam!(scene, (θ, 0, 0))
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