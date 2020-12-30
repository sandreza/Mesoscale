using JLD2, LinearAlgebra, Oceananigans
include(pwd() * "/scripts/vizinanigans.jl")
include(pwd() * "/scripts/states.jl")

include(pwd() * "/analysis_scripts/" * "post_analysis.jl") # Gradients etc. here

files = [pwd() * "/Channel_16_checkpoint_iteration6317902.jld2", 
         pwd() * "/Channel_4_checkpoint_iteration6160221.jld2",
         pwd() * "/Channel_1_checkpoint_iteration404905.jld2",
         pwd() * "/Channel_8_checkpoint_iteration3164301.jld2",
         pwd() * "/Channel_32_checkpoint_iteration1272125.jld2"
]

filename = files[end]
states, statenames = grabstates(filename)
scene = visualize(states, statenames = statenames, aspect = (1,1, 32/192), statistics = true)
display(scene)
## save interaction
seconds = 20
fps = 30
frames = round(Int, fps * seconds )
record(scene, pwd() * "/test.mp4"; framerate = fps) do io
    for i = 1:frames
        sleep(1/fps)
        recordframe!(io)
    end
end
