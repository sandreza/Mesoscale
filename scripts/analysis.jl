using JLD2, LinearAlgebra, Oceananigans, Printf, Statistics
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

new_files = [
    pwd() * "/Channel_1_checkpoint_iteration370489.jld2",
    pwd() * "/Channel_3_checkpoint_iteration1075121.jld2",
    pwd() * "/Channel_16_checkpoint_iteration1003423.jld2",
]

ridge_files = [
    pwd() * "/Ridge_8_checkpoint_iteration939858.jld2",
    pwd() * "/Ridge_8_checkpoint_iteration6078160.jld2",
    pwd() * "/Ridge_16_checkpoint_iteration1051200.jld2",
]

relaxation_files = [
    pwd() * "/Relaxation_Channel_16_checkpoint_iteration1445551.jld2"
]

abernathy_files = [
    pwd() * "/Abernathy_16_checkpoint_iteration1158859.jld2"
]
##
# http://juliaplots.org/MakieReferenceImages/gallery/index.html
filename = abernathy_files[end]
states, statenames, units = grabstates(filename)
scene = volumeslice(states, statenames = statenames, aspect = (1, 1, 32/192), 
                    statistics = true, units = units, statlabelsize = (15, 15) );
# display(scene)
## save interaction
seconds = 20
fps = 10
frames = round(Int, fps * seconds )
if record_interaction
GLMakie.record(scene, pwd() * "/ridge.mp4"; framerate = fps) do io
    for i = 1:frames
        sleep(1/fps)
        recordframe!(io)
    end
end
end
##
filename = new_files[1]
states, statenames, units1 = grabstates(filename)
title = "" * grabtitle(filename)

filename2 = new_files[end]
states2, statenames2, units2 = grabstates(filename2)
title2 = " "  * grabtitle(filename2)

scene = visualize(states, states2, statenames = statenames, statenames2 = statenames2, aspect = (1,1, 32/192), statistics = true, title = title, title2 = title2, units1 = units1, units2 = units2)
##
# instantaneous zonal 
filename = pwd() * "/Weno_20_checkpoint_iteration21030963.jld2" # new_files[1]
states, statenames, units1 = grabstates(filename)
function zonalmean(state)
    return mean(state, dims = 1)[1,:,:]
end
meanstates = zonalmean.(states)
statenames[3]
v̅ =  zonalmean(states[3])
b̅ =  zonalmean(states[5])
vb = zonalmean(states[3] .* states[5])
visualize([v̅, b̅, vb, vb - v̅ .* b̅ ], statenames = ["v̄", "b̄", "avg(vb)",  " v'b' "], title = "Zonal Average")
scene = visualize(meanstates, statenames = statenames, title = "Zonal Average ")
##
filename = new_files[end]
states, statenames, units1 = grabstates(filename)
meanstates2 = zonalmean.(states)
coarsegrainedb = coarsey(meanstates2[5], 12)
coarseb = meanstates[5]
scene = visualize([coarseb, coarsegrainedb, (coarsegrainedb-coarseb) ./ norm(coarseb)], statenames = ["b", "grained b", "difference"], title = "Coarse ")
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


