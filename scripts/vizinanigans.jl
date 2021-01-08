using Oceananigans
using GLMakie, AbstractPlotting
using ImageTransformations, Colors
using AbstractPlotting.MakieLayout
using Statistics, Printf

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
visualize(states::AbstractArray; statenames = string.(1:length(states)), quantiles = (0.1, 0.99), aspect = false, resolution = (1920, 1080), statistics = false, title = "Field = ")

# Description 
Visualize 3D states 

# Arguments
- `states`: Array{Array{Float64,3},1}. An array of arrays containing different fields

# Keyword Arguments
- `statenames`: Array{String,1}. An array of stringnames
- `aspect`: Tuple{Int64,Int64,Float64}. Determines aspect ratio of box for volumes
- `resolution`: Resolution of preliminary makie window
- `statistics`: boolean. toggle for displaying statistics 

# Return
- `scene`: Scene. A preliminary scene object for manipulation
"""
function visualize(states::AbstractArray; statenames = string.(1:length(states)), units = ["" for i in eachindex(states)], aspect = false, resolution = (1920, 1080), statistics = false, title = "Field = ", bins = 300)
    # Create scene
    scene, layout = layoutscene(resolution = resolution)
    lscene = layout[2:4, 2:4] = LScene(scene) 
    width = round(Int, resolution[1] / 4) # make menu 1/4 of preliminary resolution

    # Create choices and nodes
    stateindex = collect(1:length(states))
    statenode = Node(stateindex[1])

    colorchoices = [:balance, :thermal, :dense, :deep, :curl, :thermometer]
    colornode = Node(colorchoices[1])

    if statistics
        llscene = layout[4,1] = LAxis(scene, xlabel = @lift(statenames[$statenode] * units[$statenode]), 
                         xlabelcolor = :black, ylabel = "pdf", 
                         ylabelcolor = :black, xlabelsize = 40, ylabelsize = 40,
                         xticklabelsize = 25, yticklabelsize = 25,
                         xtickcolor = :black, ytickcolor = :black,
                         xticklabelcolor  = :black, yticklabelcolor = :black)
        layout[3, 1] = LText(scene, "Statistics", width = width, textsize = 50)
    end

    # x,y,z are for determining the aspect ratio of the box
    if (typeof(aspect) <: Tuple) & (length(aspect) == 3)
        x, y, z = aspect
    else
        x, y, z = size(states[1])
    end

    # Clim sliders
    upperclim_slider = LSlider(scene, range = range(0, 1, length = 101), startvalue = 0.99)
    upperclim_node = upperclim_slider.value
    lowerclim_slider = LSlider(scene, range = range(0, 1, length = 101), startvalue = 0.01)
    lowerclim_node = lowerclim_slider.value

    # Lift Nodes
    state = @lift(states[$statenode])
    statename = @lift(statenames[$statenode])
    clims = @lift((quantile($state[:], $lowerclim_node) , quantile($state[:], $upperclim_node)))
    cmap_rgb = @lift(to_colormap($colornode))
    titlename = @lift(title * $statename) # use padding and appropriate centering

    # Statistics
    if statistics
        histogram_node = @lift(histogram($state, bins = bins))
        xs = @lift($histogram_node[1])
        ys = @lift($histogram_node[2])
        pdf = AbstractPlotting.barplot!(llscene, xs, ys, color = :red, 
                        strokecolor = :red, 
                        strokewidth = 1)
        @lift(AbstractPlotting.xlims!(llscene, extrema($state)))
        @lift(AbstractPlotting.ylims!(llscene, extrema($histogram_node[2])))
        vlines!(llscene, @lift($clims[1]), color = :black, linewidth = width / 100)
        vlines!(llscene, @lift($clims[2]), color = :black, linewidth = width / 100)
    end

    # Volume Plot 
    volume!(lscene, 0..x, 0..y, 0..z, state, 
            camera = cam3d!, 
            colormap = cmap_rgb, 
            colorrange = clims)
    # Camera
    cam = cameracontrols(scene.children[1])
    eyeposition = Float32[2, 2, 1.3]
    lookat = Float32[0.82, 0.82, 0.1]
    # Title
    supertitle = layout[1, 2:4] = LText(scene, titlename , textsize = 50, color = :black)
    

    # Menus
    statemenu = LMenu(scene, options = zip(statenames, stateindex))
    on(statemenu.selection) do s
        statenode[] = s
        update_cam!(scene.children[1], eyeposition, lookat, Vec3f0(0, 0, 1))
    end

    colormenu = LMenu(scene, options = zip(colorchoices, colorchoices))
    on(colormenu.selection) do s
        colornode[] = s
        update_cam!(scene.children[1], eyeposition, lookat, Vec3f0(0, 0, 1))
    end
    lowerclim_string = @lift("lower clim quantile = " *  @sprintf("%0.2f", $lowerclim_node) * ", value = " * @sprintf("%0.1e", $clims[1]))
    upperclim_string = @lift("upper clim quantile = " *  @sprintf("%0.2f", $upperclim_node) * ", value = " * @sprintf("%0.1e", $clims[2]))
    # depends on makie version, vbox for old, vgrid for new
    layout[2, 1] = vgrid!(
        LText(scene, "State", width = nothing),
        statemenu,
        LText(scene, "Color", width = nothing),
        colormenu,
        LText(scene, lowerclim_string, width = nothing),
        lowerclim_slider,
        LText(scene, upperclim_string, width = nothing),
        upperclim_slider,
    )
    layout[1,1] = LText(scene, "Menu", width = width, textsize = 50)

    # Modify Axis
    axis = scene.children[1][Axis] 
    # axis[:names][:axisnames] = ("↓ Zonal [m] ", "Meriodonal [m]↓ ", "Depth [m]↓ ")
    axis[:names][:axisnames] = ("↓", "↓ ", "↓ ")
    axis[:names][:align] = ((:left, :center), (:right, :center), (:right, :center))
    # need to adjust size of ticks first and then size of axis names
    axis[:names][:textsize] = (50.0, 50.0, 50.0)
    axis[:ticks][:textsize] = (00.0, 00.0, 00.0)
    # axis[:ticks][:ranges_labels].val # current axis labels
    xticks = collect(range(-0, aspect[1], length = 2))
    yticks = collect(range(-0, aspect[2], length = 6))
    zticks = collect(range(-0, aspect[3], length = 2))
    ticks = (xticks, yticks, zticks)
    axis[:ticks][:ranges] = ticks
    xtickslabels = [@sprintf("%0.1f", (xtick)) for xtick in xticks]
    xtickslabels[end] = "1e6"
    ytickslabels = ["", "south","", "", "north", ""]
    ztickslabels = [@sprintf("%0.1f", (xtick)) for xtick in xticks]
    labels = (xtickslabels, ytickslabels, ztickslabels)
    axis[:ticks][:labels] = labels

    update_cam!(scene.children[1], eyeposition, lookat, Vec3f0(0, 0, 1))
    display(scene)
    # Change the default camera position after the fact
    # note that these change dynamically as the plot is manipulated
    update_cam!(scene.children[1], eyeposition, lookat, Vec3f0(0, 0, 1))
    return scene
end


"""
histogram(array; bins = 100)

# Description
return arrays for plotting histogram
"""
function histogram(array; bins = minimum([100, length(array)]), normalize = true)
    tmp = zeros(bins)
    down, up = extrema(array)
    down, up = down == up ? (down-1, up+1) : (down, up) # edge case
    bucket = collect(range(down, up, length = bins+1))
    normalization = normalize ? length(array) : 1
    for i in eachindex(array)
        # normalize then multiply by bins
        val = (array[i] - down) / (up - down) * bins
        ind = ceil(Int, val)
        # handle edge cases
        ind = maximum([ind, 1])
        ind = minimum([ind, bins])
        tmp[ind] += 1/normalization
    end
    return (bucket[2:end] + bucket[1:end-1]) .* 0.5, tmp
end

# 2D visualization
function visualize(states::Array{Array{S, 2},1}; statenames = string.(1:length(states)), units = ["" for i in eachindex(states)], aspect = false, resolution = (2412, 1158), title = "Zonal and Temporal Average of ", xlims = (0,1), ylims = (0,1), bins = 300) where {S}
    # Create scene
    scene, layout = layoutscene(resolution = resolution)
    lscene = layout[2:4, 2:4] = LAxis(scene, xlabel = "South to North [m]", 
                         xlabelcolor = :black, ylabel = "Depth [m]", 
                         ylabelcolor = :black, xlabelsize = 40, ylabelsize = 40,
                         xticklabelsize = 25, yticklabelsize = 25,
                         xtickcolor = :black, ytickcolor = :black,
                         xticklabelcolor  = :black, yticklabelcolor = :black,
                         titlesize = 50) 
    width = round(Int, resolution[1] / 4) # make menu 1/4 of preliminary resolution

    # Create choices and nodes
    stateindex = collect(1:length(states))
    statenode = Node(stateindex[1])

    colorchoices = [:balance, :thermal, :dense, :deep, :curl, :thermometer]
    colornode = Node(colorchoices[1])

    interpolationlabels = ["contour", "heatmap"]
    interpolationchoices = [true, false]
    interpolationnode = Node(interpolationchoices[1])

    # Statistics
    llscene = layout[4,1] = LAxis(scene, xlabel = @lift(statenames[$statenode] * " " * units[$statenode]), 
                    xlabelcolor = :black, ylabel = "pdf", 
                    ylabelcolor = :black, xlabelsize = 40, ylabelsize = 40,
                    xticklabelsize = 25, yticklabelsize = 25,
                    xtickcolor = :black, ytickcolor = :black,
                    xticklabelcolor  = :black, yticklabelcolor = :black)
    layout[3, 1] = LText(scene, "Statistics", width = width, textsize = 50)

    # Clim sliders
    upperclim_slider = LSlider(scene, range = range(0, 1, length = 101), startvalue = 0.99)
    upperclim_node = upperclim_slider.value
    lowerclim_slider = LSlider(scene, range = range(0, 1, length = 101), startvalue = 0.01)
    lowerclim_node = lowerclim_slider.value
   
    #ylims = @lift(range($lowerval, $upperval, length = $))
    # Lift Nodes
    state = @lift(states[$statenode])
    statename = @lift(statenames[$statenode])
    unit = @lift(units[$statenode])
    oclims = @lift((quantile($state[:], $lowerclim_node) , quantile($state[:], $upperclim_node)))
    cmap_rgb = @lift($oclims[1] < $oclims[2] ? to_colormap($colornode) : reverse(to_colormap($colornode)))
    clims = @lift($oclims[1] != $oclims[2] ? (minimum($oclims), maximum($oclims)) : (minimum($oclims)-1, maximum($oclims)+1))
    xlims = Array(range(xlims[1], xlims[2], length = 4)) #collect(range(xlims[1], xlims[2], length = size(states[1])[1]))
    ylims = Array(range(ylims[1], ylims[2], length = 4)) #@lift(collect(range($lowerval], $upperval, length = size($state)[2])))
    # newrange = @lift(range($lowerval, $upperval, length = 4))
    # lscene.yticks = @lift(Array($newrange))
    titlename = @lift(title * $statename) # use padding and appropriate centering
    layout[1, 2:4] = LText(scene, titlename, textsize = 50)
    # heatmap 
    heatmap1 = heatmap!(lscene, xlims, 
            ylims,
            state, interpolate = interpolationnode,
            colormap = cmap_rgb, colorrange = clims)


    # statistics
    histogram_node = @lift(histogram($state, bins = bins))
    xs = @lift($histogram_node[1])
    ys = @lift($histogram_node[2])
    pdf = AbstractPlotting.barplot!(llscene, xs, ys, color = :red, 
                    strokecolor = :red, 
                    strokewidth = 1)
    @lift(AbstractPlotting.xlims!(llscene, extrema($state)))
    @lift(AbstractPlotting.ylims!(llscene, extrema($histogram_node[2])))
    vlines!(llscene, @lift($clims[1]), color = :black, linewidth = width / 100)
    vlines!(llscene, @lift($clims[2]), color = :black, linewidth = width / 100)

    # Menus
    statemenu = LMenu(scene, options = zip(statenames, stateindex))
    on(statemenu.selection) do s
        statenode[] = s
    end

    colormenu = LMenu(scene, options = zip(colorchoices, colorchoices))
    on(colormenu.selection) do s
        colornode[] = s
    end

    interpolationmenu = LMenu(scene, options = zip(interpolationlabels, interpolationchoices))
    on(interpolationmenu.selection) do s
        interpolationnode[] = s
    heatmap1 = heatmap!(lscene, xlims, 
            ylims,
            state, interpolate = s,
            colormap = cmap_rgb, colorrange = clims)
    end

    newlabel = @lift($statename * " " * $unit)
    cbar = LColorbar(scene, heatmap1, label = newlabel)
    cbar.width = Relative(1/3)
    cbar.height = Relative(5/6)
    cbar.halign = :center
    cbar.flipaxisposition = true
    # cbar.labelpadding = -350
    cbar.labelsize = 50

    lowerclim_string = @lift("clim quantile = " *  @sprintf("%0.2f", $lowerclim_node) * ", value = " * @sprintf("%0.1e", $clims[1]))
    upperclim_string = @lift("clim quantile = " *  @sprintf("%0.2f", $upperclim_node) * ", value = " * @sprintf("%0.1e", $clims[2]))

    # depends on makie version, vbox for old, vgrid for new
    layout[2, 1] = vgrid!(
        LText(scene, "State", width = nothing),
        statemenu,
        LText(scene, "plotting options", width = width, textsize = 30, padding = (0,0, 10, 0)),
        interpolationmenu,
        LText(scene, "Color", width = nothing),
        colormenu,
        LText(scene, lowerclim_string, width = nothing),
        lowerclim_slider,
        LText(scene, upperclim_string, width = nothing),
        upperclim_slider,
    )

    layout[2:4, 5] = vgrid!(
        LText(scene, "Color Bar", width = width/2, textsize = 50, padding = (25, 0, 0, 00)),
        cbar,
    )
    layout[1,1] = LText(scene, "Menu", width = width, textsize = 50)
    display(scene)
    return scene
end