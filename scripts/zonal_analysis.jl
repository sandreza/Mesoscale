using JLD2, LinearAlgebra, Oceananigans, Printf
record_interaction = false
include(pwd() * "/scripts/vizinanigans.jl")
include(pwd() * "/scripts/states.jl")
include(pwd() * "/scripts/compare.jl")
include(pwd() * "/scripts/zonalstates.jl")

include(pwd() * "/analysis_scripts/" * "post_analysis.jl") 

files = [
    pwd() * "/Channel_1_zonal_averages.jld2",
    pwd() * "/Channel_16_zonal_averages.jld2",
]
##
file = files[2]
states, statenames, units, domain = grabzonalstates(file)
li = 16 # bottom
mval = 2 # positive
ui = length(domain[2])-mval # top  
newstates = [state[:, li:end-mval] for state in states]
xlims = (domain[1][1], domain[1][end])
ylims = (domain[2][li], domain[2][ui])
scene = visualize(newstates, statenames = statenames, xlims = xlims, ylims = ylims, units = units)

##


