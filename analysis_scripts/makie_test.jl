using AbstractPlotting
using AbstractPlotting.MakieLayout

s1 = slider(LinRange(0, 100, 100), start = 10)
s2 = slider(LinRange(-2pi, 2pi, 100))
data = lift(s2[end][:value]) do v
     map(LinRange(0, 2pi, 100)) do x
         4f0 .* Point2f0(sin(x) + (sin(x * v) .* 0.1), cos(x) + (cos(x * v) .* 0.1))
     end
end
p = GLMakie.scatter(data, markersize = s1[end][:value] )
hbox(p, vbox(s1, s2))
##
using AbstractPlotting
using AbstractPlotting.MakieLayout

scene, layout = layoutscene(resolution = (1200, 900))


ax = layout[1, 1] = LAxis(scene)
sl1 = layout[2, 1] = LSlider(scene, range = 0:0.1:100, startvalue = 3)

title = " "^30 * "Buoyancy [m/s²] at ="  #need to use sl1.value to access nodes
tmp = @lift $title * string($(sl1.value))
sl2 = layout[3, 1] = LText(scene, tmp) 
sl3 = layout[4, 1] = LSlider(scene, range = 0:0.01:10, startvalue = 7)

sl4 = layout[:, 2] = LSlider(scene, range = 0:0.01:10, horizontal = false,
    tellwidth = true, height = nothing, width = Auto())
p = GLMakie.scatter!(ax, data, markersize = sl1.value )
display(scene)
