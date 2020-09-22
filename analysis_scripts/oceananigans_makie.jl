using JLD2, Plots, GLMakie, AbstractPlotting
using ImageTransformations, Colors
filename = "/Hybrid_20_checkpoint_iteration10519002.jld2"
file = jldopen( pwd() * filename)
B = file["tracers"]["b"]["data"][2:end-1, 2:end-1, 2:end-1]
close(file)
##
T = B[:,:,23:end]
# AbstractPlotting.available_gradients() # see available colors
cmap = :linear_blue_95_50_c20_n256
cmap = :thermometer
clims = (minimum(T) * 10, maximum(T))
cmapa = RGBAf0.(to_colormap(cmap), 1);
cmapa1 = vcat(fill(RGBAf0(0,0,0,0), 10), cmapa);
cmapa2 = vcat(cmapa, fill(RGBAf0(0,0,0,0), 1));
x, y, z = size(B)
scene = volume(0..x, 0..y, 0..z, T, colorrange=clims, algorithm=:absorption, absorption=10.0f0, colormap=cmapa2, show_axis=false)
volume!(x..(2x), 0..y, 0..z, T, colorrange=clims, algorithm=:absorption, absorption=3f0, colormap=cmapa2, show_axis=false)
scene
##