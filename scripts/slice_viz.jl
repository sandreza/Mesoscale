aspect = (1, 1, 32/192)
resolution = (2230, 1042)

scene, layout = layoutscene(resolution = resolution)
volumescene = layout[2:4, 2:4] = LScene(scene)
menuwidth = round(Int, resolution[1] / 4)
layout[1,1] = LText(scene, "Menu", width = menuwidth, textsize = 50)


slice_slider = LSlider(scene, range = range(0, 1, length = 101), startvalue = 0.0)
slice_node = slice_slider.value

directionindex = [1, 2, 3]
directionnames = ["x-slice", "y-slice", "z-slice"]
directionnode = Node(directionindex[1])

statenode = Node(5)
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

volume!(volumescene, 0..aspect[1], 0..aspect[2], 0..aspect[3], state, overdraw = false, colorrange = clims, colormap = to_colormap(:balance))


v = volume!(volumescene, 0..aspect[1], 0..aspect[2], 0..aspect[3],
        planeslice, algorithm = :iso, isorange = 0.005, 
        isovalue = @lift($slice_node * aspect[$directionnode]),
        transparency = true, overdraw = false, visible = true, 
        colormap = cgrad(:viridis, alpha = 0.5), colorrange = [-1,0] )

        # to_colormap(:balance)
        
nx = @lift(size($state)[1])
ny = @lift(size($state)[2])
nz = @lift(size($state)[3])


directionmenu = LMenu(scene, options = zip(directionnames, directionindex))
on(directionmenu.selection) do s
    directionnode[] = s
    slicescene = layout[3:4, 5:6] = LAxis(scene)
    if s == 1
        sliced_state = @lift( $state[round(Int, 1 + $slice_node * ($nx-1)), 1:$ny, 1:$nz]) 
        heatmap!(slicescene, @lift([1,$ny]), @lift([1,$nz]), sliced_state, interpolate = true, colormap = to_colormap(:balance))
    elseif s == 2
        sliced_state = @lift( $state[1:$nx, round(Int, 1 + $slice_node * ($ny-1)), 1:$nz])
        heatmap!(slicescene, @lift([1,$ny]), @lift([1,$nz]), sliced_state, interpolate = true, colormap = to_colormap(:balance)) 
    else
        sliced_state = @lift( $state[1:$nx, 1:$ny, round(Int, 1 + $slice_node * ($nz-1))]) 
        heatmap!(slicescene, @lift([1,$nx]), @lift([1,$ny]), sliced_state, interpolate = true, colormap = to_colormap(:balance))
    end
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
)

display(scene)