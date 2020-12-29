using JLD2, LinearAlgebra, Oceananigans
include(pwd() * "/scripts/vizinanigans.jl")
include(pwd() * "/scripts/states.jl")

include(pwd() * "/analysis_scripts/" * "post_analysis.jl") # Gradients etc. here

files = [pwd() * "/Channel_16_checkpoint_iteration6317902.jld2", 
         pwd() * "/Channel_4_checkpoint_iteration6160221.jld2",
         pwd() * "/Channel_1_checkpoint_iteration404905.jld2"
]

filename = files[1]
states, statenames = grabstates(filename)
scene = visualize(states, statenames = statenames, aspect = (1,1, 32/192), statistics = true)
display(scene)
## save interaction
fps = 10
record(scene, pwd() * "/test.mp4"; framerate = fps) do io
    for i = 1:200
        sleep(1/fps)
        recordframe!(io)
    end
end
##
scene, layout = layoutscene()
xs, ys = histogram(states[2], bins = 300)
lscene = layout[1,1] = LScene(scene)
lscene2 = layout[2,2] = LScene(scene)
bplot= AbstractPlotting.barplot!(lscene, xs, ys, color = :red, 
                strokecolor = :red, 
                strokewidth = 1,
                )
                bplot.width = 1
                bplot.height = 1 
AbstractPlotting.barplot!(lscene2, xs, ys, color = :red, 
strokecolor = :red, 
strokewidth = 1,
xticks = 0:0.1:2)
display(scene)
