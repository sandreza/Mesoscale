
function visualize(states::AbstractArray, states2::AbstractArray; statenames = string.(1:length(states)), statenames2 = string.(1:length(states2)), aspect = false, resolution = (2880, 1080), statistics = false, title = "Field = ", title2 = "Field = ", units1 = ["" for i in eachindex(states)], units2 = ["" for i in eachindex(states)])
    # Create scene
    scene, layout = layoutscene(resolution = resolution)

    lscene = layout[2:4, 2:4] = LScene(scene) 
    lscene2 = layout[2:4, 5:7] = LScene(scene) 
    width = round(Int, resolution[1] / 6) # make menu 1/4 of preliminary resolution

    # Create choices and nodes
    stateindex = collect(1:length(states))
    statenode = Node(stateindex[1])

    colorchoices = [:balance, :thermal, :dense, :deep, :curl, :thermometer]
    colornode = Node(colorchoices[1])

    stateindex2 = collect(1:length(states2))
    statenode2 = Node(stateindex2[1])

    colorchoices2 = [:balance, :thermal, :dense, :deep, :curl, :thermometer]
    colornode2 = Node(colorchoices2[1])

    if statistics
        llscene = layout[4,1] = Axis(scene, xlabel = @lift(statenames[$statenode] * units1[$statenode]), 
                         xlabelcolor = :black, ylabel = "pdf", 
                         ylabelcolor = :black, xlabelsize = 40, ylabelsize = 40,
                         xticklabelsize = 25, yticklabelsize = 25,
                         xtickcolor = :black, ytickcolor = :black,
                         xticklabelcolor  = :black, yticklabelcolor = :black)
        layout[3, 1] = Label(scene, "Statistics", width = width, textsize = 50)

        llscene2 = layout[4,8] = Axis(scene, xlabel = @lift(statenames2[$statenode2] * units2[$statenode2]), 
                         xlabelcolor = :black, ylabel = "pdf", 
                         ylabelcolor = :black, xlabelsize = 40, ylabelsize = 40,
                         xticklabelsize = 25, yticklabelsize = 25,
                         xtickcolor = :black, ytickcolor = :black,
                         xticklabelcolor  = :black, yticklabelcolor = :black)
        layout[3, 8] = Label(scene, "Statistics", width = width, textsize = 50)
    end

    # x,y,z are for determining the aspect ratio of the box
    if (typeof(aspect) <: Tuple) & (length(aspect) == 3)
        x, y, z = aspect
    else
        x, y, z = size(states[1])
    end

    # Clim sliders
    upperclim_slider = Slider(scene, range = range(0, 1, length = 101), startvalue = 0.99)
    upperclim_node = upperclim_slider.value
    lowerclim_slider = Slider(scene, range = range(0, 1, length = 101), startvalue = 0.01)
    lowerclim_node = lowerclim_slider.value

    upperclim_slider2 = Slider(scene, range = range(0, 1, length = 101), startvalue = 0.99)
    upperclim_node2 = upperclim_slider2.value
    lowerclim_slider2 = Slider(scene, range = range(0, 1, length = 101), startvalue = 0.01)
    lowerclim_node2 = lowerclim_slider2.value

    # Lift Nodes
    state = @lift(states[$statenode])
    statename = @lift(statenames[$statenode])
    clims = @lift((quantile($state[:], $lowerclim_node) , quantile($state[:], $upperclim_node)))
    cmap_rgb = @lift(to_colormap($colornode))
    titlename = @lift(title * $statename) # use padding and appropriate centering

    state2 = @lift(states2[$statenode2])
    statename2 = @lift(statenames2[$statenode2])
    clims2 = @lift((quantile($state2[:], $lowerclim_node2) , quantile($state2[:], $upperclim_node2)))
    cmap_rgb2 = @lift(to_colormap($colornode2))
    titlename2 = @lift(title2 * $statename2) # use padding and appropriate centering

    # Statistics
    if statistics
        histogram_node = @lift(histogram($state, bins = 300))
        xs = @lift($histogram_node[1])
        ys = @lift($histogram_node[2])
        pdf = GLMakie.AbstractPlotting.barplot!(llscene, xs, ys, color = :red, 
                        strokecolor = :red, 
                        strokewidth = 1)
        @lift(GLMakie.AbstractPlotting.xlims!(llscene, extrema($state)))
        @lift(GLMakie.AbstractPlotting.ylims!(llscene, extrema($histogram_node[2])))
        vlines!(llscene, @lift($clims[1]), color = :black, linewidth = width / 100)
        vlines!(llscene, @lift($clims[2]), color = :black, linewidth = width / 100)
        # 2
        histogram_node2 = @lift(histogram($state2, bins = 300))
        xs2 = @lift($histogram_node2[1])
        ys2 = @lift($histogram_node2[2])
        pdf2 = GLMakie.AbstractPlotting.barplot!(llscene2, xs2, ys2, color = :red, 
                        strokecolor = :red, 
                        strokewidth = 1)
        @lift(GLMakie.AbstractPlotting.xlims!(llscene2, extrema($state2)))
        @lift(GLMakie.AbstractPlotting.ylims!(llscene2, extrema($histogram_node2[2])))
        vlines!(llscene2, @lift($clims2[1]), color = :black, linewidth = width / 100)
        vlines!(llscene2, @lift($clims2[2]), color = :black, linewidth = width / 100)
    end

    # Volume Plot 
    volume!(lscene, 0..x, 0..y, 0..z, state, 
            camera = cam3d!, 
            colormap = cmap_rgb, 
            colorrange = clims)
    volume!(lscene2, 0..x, 0..y, 0..z, state2, 
        camera = cam3d!, 
        colormap = cmap_rgb2, 
        colorrange = clims2)
    # Title
    supertitle = layout[1, 2:4] = Label(scene, titlename , textsize = 50, color = :black)
    supertitle2 = layout[1, 5:7] = Label(scene, titlename2 , textsize = 50, color = :black)
    

    # Menus 1
    statemenu = Menu(scene, options = zip(statenames, stateindex))
    on(statemenu.selection) do s
        statenode[] = s
    end

    colormenu = Menu(scene, options = zip(colorchoices, colorchoices))
    on(colormenu.selection) do s
        colornode[] = s
    end
    lowerclim_string = @lift("quantile = " *  @sprintf("%0.2f", $lowerclim_node) * ", value = " * @sprintf("%0.1e", $clims[1]))
    upperclim_string = @lift("quantile = " *  @sprintf("%0.2f", $upperclim_node) * ", value = " * @sprintf("%0.1e", $clims[2]))
    # depends on makie version, vbox for old, vgrid for new
    layout[2, 1] = vgrid!(
        Label(scene, "State", width = nothing),
        statemenu,
        Label(scene, "Color", width = nothing),
        colormenu,
        Label(scene, lowerclim_string, width = nothing),
        lowerclim_slider,
        Label(scene, upperclim_string, width = nothing),
        upperclim_slider,
    )
    layout[1,1] = Label(scene, "Menu", width = width, textsize = 50)

    # Menus 2
    statemenu2 = Menu(scene, options = zip(statenames2, stateindex2))
    on(statemenu2.selection) do s
        statenode2[] = s
    end

    colormenu2 = Menu(scene, options = zip(colorchoices2, colorchoices2))
    on(colormenu2.selection) do s
        colornode2[] = s
    end
    lowerclim_string2 = @lift("quantile = " *  @sprintf("%0.2f", $lowerclim_node2) * ", value = " * @sprintf("%0.1e", $clims2[1]))
    upperclim_string2 = @lift("quantile = " *  @sprintf("%0.2f", $upperclim_node2) * ", value = " * @sprintf("%0.1e", $clims2[2]))
    # depends on makie version, vbox for old, vgrid for new
    layout[2, 8] = vgrid!(
        Label(scene, "State", width = nothing),
        statemenu2,
        Label(scene, "Color", width = nothing),
        colormenu2,
        Label(scene, lowerclim_string2, width = nothing),
        lowerclim_slider2,
        Label(scene, upperclim_string2, width = nothing),
        upperclim_slider2,
    )
    layout[1,8] = Label(scene, "Menu", width = width, textsize = 50)

    ## Axis
        # Modify Axis
    axis = scene.children[1][OldAxis] 
    axis[:names][:axisnames] = ("↓ Zonal [m] ", "Meriodonal [m]↓ ", "Depth [m]↓ ")
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

        # Modify Axis
    axis = scene.children[2][OldAxis] 
    axis[:names][:axisnames] = ("↓ Zonal [m] ", "Meriodonal [m]↓ ", "Depth [m]↓ ")
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

    display(scene)
    return scene
end

function grabtitle(filename)
    resolutionnumber = parse(Int64, split(filename, "_")[end-2])
    return string(round(Int, 1000 * 16 / 192 / resolutionnumber)) * " km resolution, "
end