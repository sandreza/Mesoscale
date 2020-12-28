using JLD2
include(pwd() * "/scripts/vizinanigans.jl")
filename = pwd() * "/Channel_16_checkpoint_iteration6317902.jld2"
file = jldopen( filename )
b = file["tracers"]["b"]["data"][2:end-1, 2:end-1, 2:end-1]
u = file["velocities"]["u"]["data"][2:end-1, 2:end-1, 2:end-1]
v = file["velocities"]["v"]["data"][2:end-1, 2:end-1, 2:end-1]
w = file["velocities"]["w"]["data"][2:end-1, 2:end-1, 2:end-1]
close(file)
visualize([u,v,w, b], statenames = ["u", "v", "w", "b"])

visualize([u,v,w, b[2:end-1, 2:end-1, end-15:end-4]], statenames = ["u", "v", "w", "b"])