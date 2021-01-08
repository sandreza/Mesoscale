aspect = (1, 1, 32/192)
resolution = (2230, 1042)
bins = 300

scene, layout = layoutscene(resolution = resolution)
volumescene = layout[2:4, 2:4] = LScene(scene)
menuwidth = round(Int, resolution[1] / 4)
layout[1,1] = LText(scene, "Menu", width = menuwidth, textsize = 50)



slice_slider = LSlider(scene, range = range(0, 1, length = 101), startvalue = 0.0)
slice_node = slice_slider.value

directionindex = [1, 2, 3]
directionnames = ["x-slice", "y-slice", "z-slice"]
directionnode = Node(directionindex[1])

statenode = Node(4)
colorchoices = [:balance, :thermal, :dense, :deep, :curl, :thermometer]
colornode = Node(colorchoices[1])

state = @lift(states[$statenode])
nslider = 100
xrange = range(0.00, aspect[1], length = nslider)
yrange = range(0.00, aspect[2], length = nslider)
zrange = range(0.00, aspect[3], length = nslider)
constx = collect(reshape(xrange, (nslider,1,1)))
consty = collect(reshape(yrange, (1,nslider,1)))
constz = collect(reshape(zrange, (1,1,nslider)))
matx = zeros(nslider,nslider,nslider)
maty = zeros(nslider,nslider,nslider)
matz = zeros(nslider,nslider,nslider)
matx .= constx
maty .= consty
matz .= constz
sliceconst = [matx, maty, matz]
planeslice = @lift(sliceconst[$directionnode])

upperclim_slider = LSlider(scene, range = range(0, 1, length = 101), startvalue = 0.99)
upperclim_node = upperclim_slider.value
lowerclim_slider = LSlider(scene, range = range(0, 1, length = 101), startvalue = 0.01)
lowerclim_node = lowerclim_slider.value

clims = @lift((quantile($state[:], $lowerclim_node) , quantile($state[:], $upperclim_node)))

volume!(volumescene, 0..aspect[1], 0..aspect[2], 0..aspect[3], state, overdraw = false, colorrange = clims, colormap = @lift(to_colormap($colornode)))


v = volume!(volumescene, 0..aspect[1], 0..aspect[2], 0..aspect[3],
        planeslice, algorithm = :iso, isorange = 0.005, 
        isovalue = @lift($slice_node * aspect[$directionnode]),
        transparency = true, overdraw = false, visible = true, 
        colormap = cgrad(:viridis, alpha = 0.5), colorrange = [-1,0] )

        # to_colormap(:balance)
        
nx = @lift(size($state)[1])
ny = @lift(size($state)[2])
nz = @lift(size($state)[3])


# Slice
sliceupperclim_slider = LSlider(scene, range = range(0, 1, length = 101), startvalue = 0.99)
sliceupperclim_node = sliceupperclim_slider.value
slicelowerclim_slider = LSlider(scene, range = range(0, 1, length = 101), startvalue = 0.01)
slicelowerclim_node = slicelowerclim_slider.value
oclims = @lift((quantile($sliced_state[:], $slicelowerclim_node) , quantile($sliced_state[:], $sliceupperclim_node)))
cmap_rgb = @lift($oclims[1] < $oclims[2] ? to_colormap($colornode) : reverse(to_colormap($colornode)))
sliceclims = @lift($oclims[1] != $oclims[2] ? (minimum($oclims), maximum($oclims)) : (minimum($oclims)-1, maximum($oclims)+1))

slicexaxislabel = @lift(["y", "x", "x"][$directionnode])
sliceyaxislabel = @lift(["z", "z", "y"][$directionnode])

slicexaxis = @lift([[1,$ny], [1, $nx], [1,$nx]][$directionnode])
sliceyaxis = @lift([[1,$nz], [1, $nz], [1,$ny]][$directionnode])
slicexaxis = [0,1]
sliceyaxis = [0,1]

slicescene = layout[2:3, 5:6] = LAxis(scene, xlabel = slicexaxislabel, ylabel = sliceyaxislabel)

sliced_state1 = @lift( $state[round(Int, 1 + $slice_node * ($nx-1)), 1:$ny, 1:$nz])
sliced_state2 = @lift( $state[1:$nx, round(Int, 1 + $slice_node * ($ny-1)), 1:$nz])
sliced_state3 = @lift( $state[1:$nx, 1:$ny, round(Int, 1 + $slice_node * ($nz-1))]) 
sliced_states = @lift([$sliced_state1, $sliced_state2, $sliced_state3])
sliced_state = @lift($sliced_states[$directionnode]) 

heatmap!(slicescene, slicexaxis, sliceyaxis, sliced_state, interpolate = true, colormap = @lift(to_colormap($colornode)), colorrange = sliceclims)

sliceindex = @lift([round(Int, 1 + $slice_node * ($nx-1)), round(Int, 1 + $slice_node * ($ny-1)), round(Int, 1 + $slice_node * ($nz-1))][$directionnode])
slicestring = @lift(directionnames[$directionnode] * " at index " * @sprintf("%0.2f", $sliceindex)) 
layout[1, 5:6] = LText(scene, slicestring)


hslicescene = layout[4,5:6] = LAxis(scene, xlabel = @lift(statenames[$statenode] * " " * units[$statenode]), 
                    xlabelcolor = :black, ylabel = "pdf", 
                    ylabelcolor = :black, xlabelsize = 40, ylabelsize = 40,
                    xticklabelsize = 25, yticklabelsize = 25,
                    xtickcolor = :black, ytickcolor = :black,
                    xticklabelcolor  = :black, yticklabelcolor = :black)


# Slice Statistics

slicehistogram_node = @lift(histogram($sliced_state, bins = bins))
xs = @lift($slicehistogram_node[1])
ys = @lift($slicehistogram_node[2])
pdf = AbstractPlotting.barplot!(hslicescene, xs, ys, color = :red, 
                strokecolor = :red, 
                strokewidth = 1)

@lift(AbstractPlotting.xlims!(hslicescene, extrema($xs))) 
@lift(AbstractPlotting.ylims!(hslicescene, extrema($ys)))
vlines!(hslicescene, @lift($sliceclims[1]), color = :black, linewidth = menuwidth / 100)
vlines!(hslicescene, @lift($sliceclims[2]), color = :black, linewidth = menuwidth / 100)


interpolationmenu = LMenu(scene, options = zip(["contour", "heatmap"], [true, false]))

on(interpolationmenu.selection) do s
    interpolationnode[] = s
    # hack
    heatmap!(slicescene, slicexaxis, sliceyaxis, sliced_state, interpolate = s, colormap = @lift(to_colormap($colornode)), colorrange = sliceclims)
end

directionmenu = LMenu(scene, options = zip(directionnames, directionindex))

on(directionmenu.selection) do s
    directionnode[] = s
end

lowerclim_string = @lift("lower clim quantile = " *  @sprintf("%0.2f", $lowerclim_node) * ", value = " * @sprintf("%0.1e", $clims[1]))
upperclim_string = @lift("upper clim quantile = " *  @sprintf("%0.2f", $upperclim_node) * ", value = " * @sprintf("%0.1e", $clims[2]))

layout[2, 1] = vgrid!(
    LText(scene, lowerclim_string, width = nothing),
    lowerclim_slider,
    LText(scene, upperclim_string, width = nothing),
    upperclim_slider,
    LText(scene, "Slice Direction", width = nothing),
    directionmenu,
    LText(scene, "Slice Location", width = nothing),
    slice_slider,
    LText(scene, "Contour Plot Type", width = nothing), 
    interpolationmenu,
    LText(scene, "Slice Contour Limits", width = nothing),
    sliceupperclim_slider,
    LText(scene, "Slice Contour Limits", width = nothing), 
    slicelowerclim_slider
)

display(scene)