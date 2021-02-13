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

zstates = [vb]
zstatenames = ["vb"]
scene = visualize(zstates, statenames = zstatenames)