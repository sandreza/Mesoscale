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
    pwd() * "/Channel_24_zonal_averages.jld2",
]
##
file = files[3]
zonalstatistics = jldopen(file)
tkeys = keys(zonalstatistics["timeseries"]["t"])
states, statenames, units, domain = grabzonalstates(file)
li = 16 # bottom
mval = 2 # positive
ui = length(domain[2])-mval # top  
newstates = [state[:, li:end-mval] for state in states]
xlims = (domain[1][1], domain[1][end])
ylims = (domain[2][li], domain[2][ui])
scene = visualize(newstates, statenames = statenames, xlims = xlims, ylims = ylims, units = units)

##
seconds = 20
fps = 10
frames = round(Int, fps * seconds )
if record_interaction
record(scene, pwd() * "/zonalandtemporal.mp4"; framerate = fps) do io
    for i = 1:frames
        sleep(1/fps)
        recordframe!(io)
    end
end
end

