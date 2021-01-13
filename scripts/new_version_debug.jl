using GLMakie
# LAxis -> Axis
# LText -> Label
# LColorbar -> ColorBar
# LSlider -> Slider
# LMenu -> Menu
# AbstractPlotting. -> GLMakie.AbstractPlotting
# scene.children[1][Axis] -> scene.children[1][OldAxis] 

aspect = (1, 1, 32/192)
resolution = (1920, 1080)
statistics = true
title = "Field = "
bins = 300

scene, layout = layoutscene(resolution = resolution)
    lscene = layout[2:4, 2:4] = LScene(scene) 
    width = round(Int, resolution[1] / 4) # make menu 1/4 of preliminary resolution

    # Create choices and nodes
    stateindex = collect(1:length(states))
    statenode = Node(stateindex[1])



     colorchoices = [:balance, :thermal, :dense, :deep, :curl, :thermometer]
    colornode = Node(colorchoices[1])

    if statistics
        llscene = layout[4,1] = Axis(scene, xlabel = @lift(statenames[$statenode] * units[$statenode]), 
                         xlabelcolor = :black, ylabel = "pdf", 
                         ylabelcolor = :black, xlabelsize = 40, ylabelsize = 40,
                         xticklabelsize = 25, yticklabelsize = 25,
                         xtickcolor = :black, ytickcolor = :black,
                         xticklabelcolor  = :black, yticklabelcolor = :black)
        layout[3, 1] = Label(scene, "Statistics", width = width, textsize = 50)
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
        pdf = GLMakie.AbstractPlotting.barplot!(llscene, xs, ys, color = :red, 
                        strokecolor = :red, 
                        strokewidth = 1)
        @lift(GLMakie.AbstractPlotting.xlims!(llscene, extrema($state)))
        @lift(GLMakie.AbstractPlotting.ylims!(llscene, extrema($histogram_node[2])))
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
    supertitle = layout[1, 2:4] = Label(scene, titlename , textsize = 50, color = :black)
    

    # Menus
    statemenu = Menu(scene, options = zip(statenames, stateindex))
    on(statemenu.selection) do s
        statenode[] = s
    end

    colormenu = Menu(scene, options = zip(colorchoices, colorchoices))
    on(colormenu.selection) do s
        colornode[] = s
    end
    lowerclim_string = @lift("lower clim quantile = " *  @sprintf("%0.2f", $lowerclim_node) * ", value = " * @sprintf("%0.1e", $clims[1]))
    upperclim_string = @lift("upper clim quantile = " *  @sprintf("%0.2f", $upperclim_node) * ", value = " * @sprintf("%0.1e", $clims[2]))
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

    # Modify Axis

    axis = scene.children[1][OldAxis] 
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

    display(scene)