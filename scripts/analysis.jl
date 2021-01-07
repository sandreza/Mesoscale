using JLD2, LinearAlgebra, Oceananigans, Printf
record_interaction = false
include(pwd() * "/scripts/vizinanigans.jl")
include(pwd() * "/scripts/states.jl")
include(pwd() * "/scripts/compare.jl")

include(pwd() * "/analysis_scripts/" * "post_analysis.jl") # Gradients etc. here

files = [pwd() * "/Channel_1_checkpoint_iteration404905.jld2",
         pwd() * "/Channel_2_checkpoint_iteration728850.jld2",
         pwd() * "/Channel_3_checkpoint_iteration1087211.jld2",
         pwd() * "/Channel_4_checkpoint_iteration6160221.jld2",
         pwd() * "/Channel_8_checkpoint_iteration3164301.jld2",
         pwd() * "/Channel_16_checkpoint_iteration6317902.jld2",
         pwd() * "/Channel_24_checkpoint_iteration8612852.jld2", 
         pwd() * "/Channel_32_checkpoint_iteration1272125.jld2"
]
weak_files = [
        pwd() * "/Weak_Channel_1_checkpoint_iteration370489.jld2",
        pwd() * "/Weak_Channel_3_checkpoint_iteration1087211.jld2"
]
##
# http://juliaplots.org/MakieReferenceImages/gallery/index.html
filename = files[end-2]
states, statenames, units = grabstates(filename)
scene = visualize(states, statenames = statenames, aspect = (1,1, 32/192), statistics = true, units = units);
display(scene)
## save interaction
seconds = 20
fps = 30
frames = round(Int, fps * seconds )
if record_interaction
record(scene, pwd() * "/test.mp4"; framerate = fps) do io
    for i = 1:frames
        sleep(1/fps)
        recordframe!(io)
    end
end
end
##
filename = files[4]
states, statenames, units1 = grabstates(filename)
title = "" * grabtitle(filename)

filename2 = weak_files[2]
states2, statenames2, units2 = grabstates(filename2)
title2 = "Weak "  * grabtitle(filename2)

scene = visualize(states, states2, statenames = statenames, statenames2 = statenames2, aspect = (1,1, 32/192), statistics = true, title = title, title2 = title2, units1 = units1, units2 = units2)
##
seconds = 20
fps = 30
frames = round(Int, fps * seconds )
if record_interaction
record(scene, pwd() * "/test.mp4"; framerate = fps) do io
    for i = 1:frames
        sleep(1/fps)
        recordframe!(io)
    end
end
end


