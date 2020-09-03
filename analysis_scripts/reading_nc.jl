using NetCDF, Plots, Printf

filename = pwd() * sub_directory * "_middepth.nc"
ncinfo(filename)
x = Array(NetCDF.open(filename, "xC"))
y = Array(NetCDF.open(filename, "yC"))
z = Array(NetCDF.open(filename, "zC"))
t = Array(NetCDF.open(filename, "time"))
sim_day = t ./ 86400
z = (z[2:end] + z[1:end-1]) / 2
b = NetCDF.open(filename, "b")
u = NetCDF.open(filename, "u")
v = NetCDF.open(filename, "v")
w = NetCDF.open(filename, "w")
##
sti = 150 # starting time index
cmax = maximum(Array(b[2:end-1,2:end-1, 1, sti:end]))
cmin = minimum(Array(b[:,:, 1, sti:end]))
clims = (cmin, cmax)
anim = @animate for i in sti:1:length(b[1, 1, 1, :])
    b_array = Array(b[:, :, 1, i])
    # p1 = contourf(y, z, b_array, fill = true, linewidth = 0, color = :ocean, clim = (-0.0013871555898403098, -3.3776441214941526e-6))
    day_label = @sprintf("%.2f ", sim_day[i])
    p1 = contourf(x, y, b_array', 
    color = :ocean, title = "Buoyancy at Middepth at day " * day_label ,
     xlabel = "Meridional [m]", ylabel = "Depth [m]"
     , clims = clims)
end
gif(anim, filename *".gif", fps = 15)
##
filename = pwd() * sub_directory * "_zonal_average.nc"
ncinfo(filename)
y = Array(NetCDF.open(filename, "yC"))
z = Array(NetCDF.open(filename, "zC"))
b = NetCDF.open(filename, "Bz")
u = NetCDF.open(filename, "Uz")
t = Array(NetCDF.open(filename, "time"))
sim_day = t ./ 86400
cmax = maximum(Array(b[2:end-1, 2:end-1, 1:1:end]))
cmin = minimum(Array(b[2:end-1, 2:end-1, 1:1:end]))
clims = (cmin, cmax)
anim = @animate for i in 1:1:length(b[ 1, 1, :])
    b_array = Array(b[:, :, i])
    day_label = @sprintf("%.2f ", sim_day[i])
    p1 = contourf(y, z, b_array', 
        color = :ocean, title = "Zonal Average Buoyancy at day " * day_label ,
        xlabel = "Meridional [m]", ylabel = "Depth [m]"
        , clims = clims)
end
gif(anim, filename * ".gif", fps = 15)
##
i = length(b[ 1, 1, :])
u_array = Array(u[:, :, i])
day_label = @sprintf("%.2f ", sim_day[i])
p1 = contourf(y, z, u_array', 
color = :ocean, title = "Zonal Average Zonal Velocity at day " * day_label ,
        xlabel = "Meridional [m]", ylabel = "Depth [m]")

##
i = length(b[1, 1, :])
b_array = Array(b[:, :, i])
cmax = maximum(b_array)
cmin = minimum(b_array)
clims = (cmin, cmax)
day_label = @sprintf("%.2f ", sim_day[i])
p1 = contourf(y, z, b_array', 
        color = :ocean, title = "Zonal Average Buoyancy at day " * day_label ,
        xlabel = "Meridional [m]", ylabel = "Depth [m]"
        , clims = clims)
plot(p1)
##
filename = pwd() * sub_directory * "_surface.nc"
ncinfo(filename)
y = Array(NetCDF.open(filename, "yC"))
x = Array(NetCDF.open(filename, "xC"))
b = NetCDF.open(filename, "b")
t = Array(NetCDF.open(filename, "time"))
sim_day = t ./ 86400
start_time = 1
end_time = 220
cmax = maximum(Array(b[2:end-1, 2:end-1, 1, start_time:end_time]))
cmin = minimum(Array(b[2:end-1, 2:end-1, 1, start_time:end_time]))
clims = (cmin, cmax)
anim = @animate for i in start_time:end_time
    b_array = Array(b[:, :,1, i])
    day_label = @sprintf("%.2f ", sim_day[i])
    p1 = contourf(x, y, b_array', 
        color = :thermometer, title = "Surface Buoyancy at " * day_label ,
        xlabel = "Zonal [m]", ylabel = "Meridional [m]"
        , clims = clims, linewidth = 0)
end
gif(anim, filename * ".gif", fps = 15)
##
filename = pwd() *  "/Windstress_Convection_4_meridional.nc"
ncinfo(filename)
x = Array(NetCDF.open(filename, "yC"))
y = Array(NetCDF.open(filename, "yC"))
z = Array(NetCDF.open(filename, "zC"))
t = Array(NetCDF.open(filename, "time"))
sim_day = t ./ 86400
# z = (z[2:end] + z[1:end-1]) / 2
b = NetCDF.open(filename, "b")
u = NetCDF.open(filename, "u")
v = NetCDF.open(filename, "v")
w = NetCDF.open(filename, "w")
Δz = z[2:end] - z[1:end-1]
Δx = x[2:end] - x[1:end-1]
Δy = y[2:end] - y[1:end-1]
a = abs.(w[1, 1, 2:end-1, 1] ./ (Δz))
a2 = abs.(w[1, 1:end-1, 16, 1] ./ (Δx))
b = abs.(u[1, 1:end-1, 1, 1] ./ (Δy))
c = abs.(v[1, 2:end-1, 1, 1] ./ (Δy))

##
filename = pwd() *  "/Windstress_Convection_4_checkpoint_iteration4851612.jld2"
ncinfo(filename)
