include("local_diffusivity_yz.jl")
using GLMakie 

K = diffusivities[1]
interpolation_flag = true
fig = Figure()
ax11 = fig[1,1] = Axis(fig, title = "Kʸʸ")
field = K[1,1,:,:];
clims = quantile.(Ref(field[:]), [0.1, 0.9])
hm11 = heatmap!(ax11, field, colorrange = clims, colormap = :thermal, interpolate = interpolation_flag)

Colorbar(fig[1,2], hm11, label = " ",
topspinevisible = true, 
bottomspinevisible = true, 
leftspinevisible = true,
rightspinevisible = true)

ax12 = fig[1,3] = Axis(fig, title = "Kʸᶻ")
field = K[1,2,:,:];
clims = quantile.(Ref(field[:]), [0.1, 0.9])
hm12 = heatmap!(ax12, field, colorrange = clims, colormap = :thermal, interpolate = interpolation_flag)

Colorbar(fig[1,4], hm12, label = " ",
topspinevisible = true, 
bottomspinevisible = true, 
leftspinevisible = true,
rightspinevisible = true)

ax21 = fig[2,1] = Axis(fig, title = "Kᶻʸ")
field = K[2,1,:,:];
clims = quantile.(Ref(field[:]), [0.1, 0.9])
hm21 = heatmap!(ax21, field, colorrange = clims, colormap = :thermal, interpolate = interpolation_flag)

Colorbar(fig[2,2], hm21, label = " ",
topspinevisible = true, 
bottomspinevisible = true, 
leftspinevisible = true,
rightspinevisible = true)

ax22 = fig[2,3] = Axis(fig, title = "Kᶻᶻ")
field = K[2,2,:,:];
clims = quantile.(Ref(field[:]), [0.1, 0.9])
hm22 = heatmap!(ax22, field, colorrange = clims, colormap = :thermal, interpolate = interpolation_flag)

Colorbar(fig[2,4], hm22, label = " ",
topspinevisible = true, 
bottomspinevisible = true, 
leftspinevisible = true,
rightspinevisible = true)

