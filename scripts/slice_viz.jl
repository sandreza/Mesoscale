aspect = (1, 1, 32/192)
resolution = (1502, 1018)

scene, layout = layoutscene(resolution = resolution)
volumescene = layout[2:4, 2:4] = LScene(scene)
# slicescene = layout[1, 2] = LAxis(scene) 
menuwidth = round(Int, resolution[1] / 4)
layout[1,1] = LText(scene, "Menu", width = menuwidth, textsize = 50)


slice_slider = LSlider(scene, range = range(0, 1, length = 101), startvalue = 0.0)
slice_node = slice_slider.value

directionindex = [1, 2, 3]
directionnames = ["x-slice", "y-slice", "z-slice"]
directionnode = Node(directionindex[1])

state = states[1]
nslider = 100
xrange = range(0.01, aspect[1], length = nslider)
yrange = range(0.01, aspect[2], length = nslider)
zrange = range(0.01, aspect[3], length = nslider)
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

volume!(volumescene, 0..aspect[1], 0..aspect[2], 0..aspect[3], state, overdraw = false)


v = volume!(volumescene, 0..aspect[1], 0..aspect[2], 0..aspect[3],
        planeslice, algorithm = :iso, isorange = 0.005, 
        isovalue = @lift($slice_node * aspect[$directionnode]),
        transparency = true, overdraw = false, visible = true, 
        colormap = cgrad(:viridis, alpha = 0.5), colorrange = [-1.0, 0.0] )

        # to_colormap(:balance)
# heatmap!(slicescene, state[:,:,1])

directionmenu = LMenu(scene, options = zip(directionnames, directionindex))
on(directionmenu.selection) do s
    directionnode[] = s
end

layout[2, 1] = vgrid!(
    LText(scene, "Slice Direction", width = nothing),
    directionmenu,
    LText(scene, "Slice Location", width = nothing),
    slice_slider,
)

display(scene)