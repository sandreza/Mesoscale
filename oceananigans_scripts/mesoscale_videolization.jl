using GLMakie
using JLD2
using Statistics

videofiles = filter(x -> contains(x, "video"), readdir(pwd()))

# need to sort the files
iteration_number = zeros(Int64, length(videofiles))
for (i, videofile) in enumerate(videofiles)
    iteration_number[i] = parse(Int64, videofile[21:end-5])
end

sorted_videofiles = videofiles[sortperm(iteration_number)]
jl_file_base = jldopen(sorted_videofiles[1])
time_node = Node(1)
jl_file = @lift(jldopen(sorted_videofiles[$time_node]))
ghost = 3
# u = jl_file["velocities"]["u"]["data"][ghost+1:end-ghost, ghost+1:end-ghost, ghost+1:end-ghost]
v = @lift($jl_file["velocities"]["v"]["data"][ghost+1:end-ghost, ghost+1:end-ghost, ghost+1:end-ghost])
b = @lift($jl_file["tracers"]["b"]["data"][ghost+1:end-ghost, ghost+1:end-ghost, ghost+1:end-ghost])


# state = @lift($v .- sum($v, dims = 3) / 32)
state = b

clims = 0.25 # quantile(v[:], 0.99)
clims = (-clims, clims)
clims = (1.711131836628949e-6, 0.015225687990137533)
title_string = "Channel"
fig = Figure(resolution = (1520, 980))
ax = LScene(fig[1, 1], scenekw = (camera = cam3d!, show_axis = true))
ax_text = Label(fig[1, 1], title_string,
    textsize = 30, color = (:black, 0.85))

cmap = :afmhot # :balance # :Blues_9
cmapa = reverse(RGBAf0.(to_colormap(cmap), 1));
# cmap = vcat(cmapa[1:15], fill(RGBAf0(0, 0, 0, 0), 10), cmapa[25:end])
# cmap = vcat(cmapa[1:4:20], fill(RGBAf0(0, 0, 0, 0), 20), cmapa[25:end])
# , cmapa[8:12], fill(RGBAf0(0, 0, 0, 0), 4), cmapa[16:20], fill(RGBAf0(0, 0, 0, 0), 4), cmapa[20:24]
cmap = vcat(cmapa[1:4])
for i in 1:5
    global cmap = vcat(cmap, fill(RGBAf0(0, 0, 0, 0), 16), cmapa[1+4*i:8*i])
end

v1 = volume!(ax, 0 .. 10, 0 .. 20, 0 .. 5, state,
    colorrange = clims, algorithm = :absorption, absorption = 20.0f0,
    colormap = cmap)
axis = ax.scene[OldAxis]
axis[:names, :axisnames] = ("longitude [ᵒ]", "latitude [ᵒ]", "Depth [km]")
tstyle = axis[:names] #  get the nested attributes and work directly with them

tstyle[:textsize] = 15
tstyle[:textcolor] = (:black, :black, :black)
tstyle[:font] = "helvetica"
tstyle[:gap] = 10
axis[:ticks][:textcolor] = :black
axis[:ticks][:textsize] = 10
cbar1 = Colorbar(fig[1, 2], v1, label = " b [m/s²]", width = 25, ticklabelsize = 30,
    labelsize = 30, ticksize = 25, tickalign = 1, height = Relative(3 / 4)
)

axis[:ticks][:ranges] = ([0.0, 2.5, 5.0, 7.5, 10.0], [0.0, 5.0, 10.0, 15.0, 20.0], [0, 5 / 3, 10 / 3, 15 / 3])
axis[:ticks][:labels] = (["110E", "112E", "115E", "117E", "120E"], ["60S", "55S", "50S", "45S", "40S"], ["3", "2", "1", "0"])

fig[2:10, 1:10] = ax
fig[3:8, 11] = cbar1
fig[1, 5:6] = ax_text

##


update!(fig.scene)

display(fig)
zoom!(ax.scene, 0.8)
zoom!(ax.scene, 0.8)

#=
for i in 1:100
    time_node[] = i
    sleep(0.05)
end
=#

# bring up the figure and then run this piece of code
framerate = 30
record(fig, "buoyancy_mesoscale.mp4", 1:361, framerate = framerate) do t
    if t == 1
        zoom!(ax.scene, 0.8)
        zoom!(ax.scene, 0.8)
    end
    time_node[] = t
    println("at time t= ", t)
end