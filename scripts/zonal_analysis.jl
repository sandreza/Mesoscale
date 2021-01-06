using JLD2, LinearAlgebra, Oceananigans, Printf
record_interaction = false
include(pwd() * "/scripts/vizinanigans.jl")
include(pwd() * "/scripts/vizinanigans_2D.jl")
include(pwd() * "/scripts/states.jl")
include(pwd() * "/scripts/compare.jl")
include(pwd() * "/scripts/zonalstates.jl")

include(pwd() * "/analysis_scripts/" * "post_analysis.jl") 

files = [
    pwd() * "/Channel_1_zonal_averages.jld2",
    pwd() * "/Channel_16_zonal_averages.jld2",
]

file = files[2]
states, statenames, units = grabzonalstates(file)
scene = visualize2D(states, statenames = statenames, xlims = (0, 1e6), ylims = (-3000, 0), units = units)


