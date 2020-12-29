using Oceananigans
using GLMakie, AbstractPlotting
using ImageTransformations, Colors
using AbstractPlotting.MakieLayout
using Statistics

"""
visualize(model::Oceananigans.AbstractModel)

# Description
Create visualize tracers and velocity fields from an Oceananigans model 

# Return
- `scene`: Scene. A makie scene object
"""
function visualize(model::Oceananigans.AbstractModel)
    vstates = [Array(interior(model.velocities[velocity])) for velocity in keys(model.velocities)]
    vstatenames = [string(velocity) for velocity in keys(model.velocities)]
    tstates = [Array(interior(model.tracers[tracer])) for tracer in keys(model.tracers)]
    tstatenames = [string(tracer) for tracer in keys(model.tracers)]
    states = vcat(vstates, tstates)
    statenames = vcat(vstatenames, tstatenames)
    scene = visualize(states, statenames = statenames)
    return scene
end

"""
visualize(states::AbstractArray; statenames = string.(1:length(states)), quantiles = (0.1, 0.99), aspect = false, resolution = (1920, 1080))

# Description 
Visualize 3D states 

# Arguments
- `states`: Array{Array{Float64,3},1}. An array of arrays containing different fields

# Keyword Arguments
- `statenames`: Array{String,1}. An array of stringnames
- `quantiles`: Tuple{Number, Number}. Tuple for determining upper and lower bounds dynamically from 
- `aspect`: Tuple{Int64,Int64,Float64}. Determines aspect ratio of box for volumes
- `resolution`: Resolution of preliminary makie window

# Return
- `scene`: Scene. A preliminary scene object for manipulation
"""
function visualize(states::AbstractArray; statenames = string.(1:length(states)), quantiles = (0.1, 0.99), aspect = false, resolution = (1920, 1080))
    # Create choices and nodes
    stateindex = collect(1:length(states))
    statenode = Node(stateindex[1])

    colorchoices = [:balance, :thermal, :dense, :deep, :curl, :thermometer]
    colornode = Node(colorchoices[1])

    # x,y,z are for determining the aspect ratio of the box
    if (typeof(aspect) <: Tuple) & (length(aspect) == 3)
        x, y, z = aspect
    else
        x, y, z = size(states[1])
    end

    # Lift Nodes
    state = @lift(states[$statenode])
    clims = @lift((quantile(states[$statenode][:], quantiles[1]) , quantile(states[$statenode][:], quantiles[2]))) # lower bound not working
    cmap_rgb = @lift(to_colormap($colornode))
    titlename = @lift(" Field =" * statenames[$statenode]) # use padding and appropriate centering

    # Create scene
    scene, layout = layoutscene(resolution = resolution)
    # Volume Plot (needs to come first)
    lscene = layout[2:4, 2:4] = LScene(scene) 
    volume!(lscene, 0..x, 0..y, 0..z, state, 
            camera = cam3d!, 
            colormap = cmap_rgb, 
            colorrange = clims)
    # Title
    supertitle = layout[1,2:4] = LText(scene, titlename , textsize = 50, color = :black)
    # Menus
    statemenu = LMenu(scene, options = zip(statenames, stateindex))
    on(statemenu.selection) do s
        statenode[] = s
    end

    colormenu = LMenu(scene, options = zip(colorchoices, colorchoices))
    on(colormenu.selection) do s
        colornode[] = s
    end

    layout[1, 1] = vgrid!(
        LText(scene, "State", width = nothing),
        statemenu,
        LText(scene, "Color", width = nothing),
        colormenu,
    )
    display(scene)
    return scene
end
