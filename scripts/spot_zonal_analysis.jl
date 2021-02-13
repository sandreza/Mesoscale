using JLD2, LinearAlgebra, Oceananigans, Printf
record_interaction = false
include(pwd() * "/scripts/vizinanigans.jl")
include(pwd() * "/scripts/states.jl")
include(pwd() * "/scripts/compare.jl")
include(pwd() * "/scripts/zonalstates.jl")
include(pwd() * "/analysis_scripts/" * "post_analysis.jl") 

file = pwd() * "/Abernathy_16_zonal_averages.jld2"
zonalstatistics = jldopen(file)
tkeys = keys(zonalstatistics["timeseries"]["t"])

i = 10

label = "vb"
vb = zonalstatistics["timeseries"][label][tkeys[i]][1,:,:]

##
filename = pwd() * "/Abernathy_16_checkpoint_iteration1158859.jld2"
states, statenames, units = grabstates(filename)
v = states[3]
b = states[5]
vb = mean(v .* b, dims = 1)[1,:,:]
vpbp = vb - mean(v, dims=1)[1,:,:] .* mean(b, dims=1)[1,:,:]
zstates = [vb, vpbp]
zstatenames = ["vb", "v'b'"]
scene = visualize(zstates, statenames = zstatenames)