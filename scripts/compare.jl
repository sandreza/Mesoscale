files = [pwd() * "/Channel_16_checkpoint_iteration6317902.jld2", 
         pwd() * "/Channel_4_checkpoint_iteration6160221.jld2",
         pwd() * "/Channel_1_checkpoint_iteration404905.jld2"
]
statenames = ["u"]
stateindex = [1]
x = y = 1
z = 32/192
cmap_rgb = :balance
colorchoices = [:balance]
state = u
# x, y, z = size(state) # basically only determines the aspect ratio
minfield, maxfield = extrema(state)
# Lift Nodes
state = (state .- minfield) ./ (maxfield - minfield)
(lowq, upperq) = (quantile(state[:], 0.1), quantile(state[:], 0.99))
clims = (lowq, upperq) 
titlename = " Field "# use padding and appropriate centering
# Create scene
scene, layout = layoutscene(resolution = (1920, 860))
# Volume Plot (needs to come first)
lscene = layout[2, 2:3] = LScene(scene) 
volume!(lscene, 0..x, 0..y, 0..z, state, camera = cam3d!, colormap = cmap_rgb, 
        colorrange = clims)
lscene2 = layout[2, 4:5] = LScene(scene) 
volume!(lscene2, 0..x, 0..y, 0..z, state, camera = cam3d!, colormap = cmap_rgb, 
        colorrange = (1,0))

lscene3 = layout[4, 2:3] = LScene(scene) 
volume!(lscene3, 0..x, 0..y, 0..z, state, camera = cam3d!, colormap = cmap_rgb, 
        colorrange = (1,0))

# Title
supertitle = layout[1,2:3] = LText(scene, titlename , textsize = 50, color = :black)
supertitle2 = layout[1,4:5] = LText(scene, titlename , textsize = 50, color = :black)

supertitle3 = layout[3,4:5] = LText(scene, titlename , textsize = 50, color = :black)

statemenu = LMenu(scene, options = zip(statenames, stateindex))
on(statemenu.selection) do s
    statenode[] = s
end

colormenu = LMenu(scene, options = zip(colorchoices, colorchoices))
on(colormenu.selection) do s
    colornode[] = s
end

statemenu2 = LMenu(scene, options = zip(statenames, stateindex))
on(statemenu2.selection) do s
    statenode[] = s
end

colormenu2 = LMenu(scene, options = zip(colorchoices, colorchoices))
on(colormenu2.selection) do s
    colornode[] = s
end

slidermenu = LSlider(scene, range = [0.1, 0.2, 0.3, 1.0], startvalue = 1)
slidernode = slidermenu.value

layout[1, 1] = vgrid!(
    LText(scene, "State", width = nothing),
    statemenu,
    LText(scene, "Color", width = nothing),
    colormenu,
        LText(scene, "Color", width = nothing),
    colormenu,
        LText(scene, "Color", width = nothing),
    colormenu,
    LText(scene, @lift("Upper Clim = " * string($slidernode)), width = nothing),
    slidermenu,
)

#sl1 = layout[2, 1] = LSlider(scene, range = 1:10, startvalue = 1)
#obs = sl1.value # make index a "node" to allow for movie making

layout[3, 6] = vgrid!(
    LText(scene, "State", width = nothing),
    statemenu2,
    LText(scene, "Color", width = nothing),
    colormenu2,
)
display(scene)