using JLD2
using GLMakie

jl_file    = jldopen("fluxernathy_averages.jld2", "r+")
# jl_file_nh = jldopen("fluxernathy_nh_averages.jld2", "r+")
jl_file_nh = jldopen("f_nh.jld2", "r+")
jl_file_nh_ns = jldopen("nh_ns_averages.jld2", "r+")

s_keys = keys(jl_file["timeseries"])
t_keys = keys(jl_file["timeseries"][s_keys[1]])

s_index = 2
t_index = 13 # 10 year averages, lets look at year 100ish

state_name = s_keys[s_index]
state_name = "u"

field = jl_file["timeseries"][state_name][t_keys[t_index]][1,:,:]
field_nh = jl_file_nh["timeseries"][state_name][t_keys[t_index]][1,:,:]
field_nh_ns = jl_file_nh_ns["timeseries"][state_name][t_keys[t_index]][1,:,:]

state_name_2 = "b"
α  = 2e-4     # [K⁻¹] thermal expansion coefficient 
g  = 9.8061   # [m/s²] gravitational constant
# misleading, converting to temperature
field_b = jl_file["timeseries"][state_name_2][t_keys[t_index]][1,:,:] ./ (α * g)
field_nh_b = jl_file_nh["timeseries"][state_name_2][t_keys[t_index]][1,:,:] ./ (α * g)
field_nh_ns_b = jl_file_nh_ns["timeseries"][state_name_2][t_keys[t_index]][1,:,:] ./ (α * g)
b_contours = collect(0.25:0.5:7.25)
# b_contours = collect(range(extrema(field_b)..., length = 10))

function get_grid(field, jl_file; ghost = 3)
    yC = jl_file["grid"]["yᵃᶜᵃ"][ghost+1:end-ghost]
    zC = jl_file["grid"]["zᵃᵃᶜ"][ghost+1:end-ghost]
    yF = jl_file["grid"]["yᵃᶠᵃ"][ghost+1:end-ghost]
    zF = jl_file["grid"]["zᵃᵃᶠ"][ghost+1:end-ghost]

    y = size(field)[1] == size(yC)[1] ? yC : yF
    z = size(field)[2] == size(zC)[1] ? zC : zF

    return y,z
end

y,z = get_grid(field, jl_file)

colorrange = extrema(field)

fig = Figure(resolution = (2024, 929))
ax_h = fig[1,1] = Axis(fig)
ax_h.title = "hydrostatic (stretched)" * state_name
ax_h.titlesize = 40
hm = heatmap!(ax_h, y, z, field, colorrange = colorrange, colormap = :thermometer, interpolate = true)
contour!(ax_h, y, z, field_b, levels = b_contours, color = :black, linewidth = 3.0,)

ax_h.ylabel = "Depth [m]"
ax_h.xlabel = "South to North [m]"
ax_h.xlabelsize = 25
ax_h.ylabelsize = 25 
ax_h.limits = (extrema(y)..., extrema(z)...)


ax_nh = fig[1,2] = Axis(fig)
ax_nh.title = "nonhydrostatic (stretched)" * state_name
ax_nh.titlesize = 40
heatmap!(ax_nh, y, z, field_nh, colorrange = colorrange, colormap = :thermometer, interpolate = true)
contour!(ax_nh, y, z, field_nh_b, levels = b_contours, color = :black, linewidth = 3.0,)
#=
ax_nh_ns = fig[1,3] = Axis(fig)
ax_nh_ns.title = "nonhydrostatic " * state_name
ax_nh_ns.titlesize = 40
zC_ns = jl_file_nh_ns["grid"]["zC"][3+1:end-3]
heatmap!(ax_nh_ns, y, zC_ns, field_nh_ns, colorrange = colorrange, colormap = :thermometer, interpolate = true)
contour!(ax_nh_ns, y, zC_ns, field_nh_ns_b, levels = b_contours, color = :black, linewidth = 3.0,)

Colorbar(fig[1,4], hm, label = state_name) # , ticks = contour_levels)




for tmpax in [ax_h, ax_nh, ax_nh_ns]
    tmpax.ylabel = "Depth [m]"
    tmpax.xlabel = "South to North [m]"
    tmpax.xlabelsize = 25
    tmpax.ylabelsize = 25 
    tmpax.limits = (extrema(y)..., extrema(z)...)
end
hideydecorations!(ax_nh)
hideydecorations!(ax_nh_ns)
=#
# ax1 = fig[jj,ii] = Axis(fig)

display(fig)

#=
ax.xticks = ([-80, -60, -30, 0, 30, 60, 80], ["80S", "60S", "30S", "0", "30N", "60N", "80N"])
pressure_levels = [1000, 850, 700, 550, 400, 250, 100, 10]
ax.yticks = (pressure_levels .* 1e2, string.(pressure_levels))
ax.yreversed = true
=#
