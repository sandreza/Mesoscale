using Oceananigans
using GLMakie, AbstractPlotting
using ImageTransformations, Colors
using AbstractPlotting.MakieLayout
using Statistics

function visualize(model::Oceananigans.AbstractModel)
    vstates = [Array(interior(model.velocities[velocity])) for velocity in keys(model.velocities)]
    vstatenames = [string(velocity) for velocity in keys(model.velocities)]
    tstates = [Array(interior(model.tracers[tracer])) for tracer in keys(model.tracers)]
    tstatenames = [string(tracer) for tracer in keys(model.tracers)]
    states = vcat(vstates, tstates)
    statenames = vcat(vstatenames, tstatenames)
    visualize(states, statenames = statenames)
    return nothing
end

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
    titlename = @lift(" "^10 * " Field =" * statenames[$statenode] * " "^10) # use padding and appropriate centering

    # Create scene
    scene, layout = layoutscene(resolution = resolution)
    # Volume Plot (needs to come first)
    lscene = layout[1:4, 2:4] = LScene(scene) 
    volume!(lscene, 0..x, 0..y, 0..z, state, 
            camera = cam3d!, 
            colormap = cmap_rgb, 
            colorrange = clims)
    # Title
    supertitle = layout[1,2] = LText(scene, titlename , textsize = 50, color = :black)
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
