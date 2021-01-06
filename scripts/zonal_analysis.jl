using JLD2, LinearAlgebra, Oceananigans, Printf
record_interaction = false
include(pwd() * "/scripts/vizinanigans.jl")
include(pwd() * "/scripts/vizinanigans_2D.jl")
include(pwd() * "/scripts/states.jl")
include(pwd() * "/scripts/compare.jl")

include(pwd() * "/analysis_scripts/" * "post_analysis.jl") 

files = [
    pwd() * "/Channel_1_zonal_averages.jld2",
    pwd() * "/Channel_16_zonal_averages.jld2",
]

file = files[2]
ghost = 3 # since using WENO
zonalstatistics = jldopen(file)
tkeys = keys(zonalstatistics["timeseries"]["t"])
t = [zonalstatistics["timeseries"]["t"][tkey] for tkey in tkeys]
yC = zonalstatistics["grid"]["yC"][ghost+1:end-ghost]
zC = zonalstatistics["grid"]["zC"][ghost+1:end-ghost]
yF = zonalstatistics["grid"]["yC"][ghost+1:end-ghost]
zF = zonalstatistics["grid"]["zC"][ghost+1:end-ghost]
fields = [:u, :v, :w, :b, :vb, :vv, :wb, :ww]
for field in fields
    label = string(field)
    @eval $field = sum([zonalstatistics["timeseries"][$label][tkey][1,:,:] for tkey in tkeys])
end
close(zonalstatistics)
##
vpbp =  average2(average1(vb) - average1(v) .* b)
bz = Δ2(b) ./ reshape((zC[2:end] - zC[1:end-1]), (1, length(zC)-1))
#bz = sum(bz, dims = 1) ./ size(bz)[1]
zCA = (zC[2:end] + zC[1:end-1])/2

fieldname = "v'b'"
xlims = (0,1e6)
ylims = (-3e3, 0)
scene, layout = layoutscene(resolution = (1920, 680))
lscene = layout[1,1] = LAxis(scene, xlabel = "South to North [m]", 
                         xlabelcolor = :black, ylabel = "Depth [m]", 
                         ylabelcolor = :black, xlabelsize = 40, ylabelsize = 40,
                         xticklabelsize = 25, yticklabelsize = 25,
                         xtickcolor = :black, ytickcolor = :black,
                         xticklabelcolor  = :black, yticklabelcolor = :black,
                         title = "Zonal and Temporal Average " * fieldname,
                         titlesize = 50)
# heatmap!(lscene, yC, zCA, vpbp / bz, interpolate = true)
xnode = Node(xlims[2])
# lscene.xticks = collect(range(xlims[1], xlims[2], length = 3))
newrange = @lift(range($lowerval, $upperval, length = 2))
lscene.xticks = ([0,1], ["0", "1"]) 
lscene.limits = GeometryBasics.HyperRectangle{2,Float32}(Float32[0.0, -3000.0], Float32[1.0f6, 3000.0])
tmp = heatmap!(lscene, 
            range(xlims[1], xlims[2], length = 2), 
            range(ylims[1], ylims[2], length = 2), 
            bz, interpolate = true, colormap=:balance)
display(scene)
##
states = [eval(field) for field in fields]
statenames = [string(field) for field in fields]
push!(states, vpbp)
push!(statenames, "v'b'")
scene = visualize2D(states, statenames = statenames, xlims = (0, 1e6), ylims = (-3000, 0))

##
upperval = Node(2)
lowerval = Node(1)

