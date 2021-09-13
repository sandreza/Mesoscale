using JLD2
using GLMakie

jl_file    = jldopen("fluxernathy_averages.jld2", "r+")
jl_file_nh = jldopen("fluxernathy_nh_averages.jld2", "r+")

s_keys = keys(jl_file["timeseries"])
t_keys = keys(jl_file["timeseries"][s_keys[1]])

s_index = 2
t_index = 30 # 10 year averages, lets look at year 100ish

state_name = s_keys[s_index]
state_name = "u"

field = jl_file["timeseries"][state_name][t_keys[t_index]][1,:,:]
field_nh = jl_file_nh["timeseries"][state_name][t_keys[t_index]][1,:,:]

state_name = "b"

field_b = jl_file["timeseries"][state_name][t_keys[t_index]][1,:,:]
field_nh_b = jl_file_nh["timeseries"][state_name][t_keys[t_index]][1,:,:]

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

fig = Figure(resolution = (1186, 596))
ax_h = fig[1,1] = Axis(fig)
hm = heatmap!(ax_h, y, z, field, colorrange = colorrange, colormap = :balance, interpolate = true)
contour!(ax_h, y, z, field_b, levels = 20, color = :black, linewidth = 3.0,)

ax_nh = fig[1,2] = Axis(fig)
heatmap!(ax_nh, y, z, field_nh, colorrange = colorrange, colormap = :balance, interpolate = true)
contour!(ax_nh, y, z, field_nh_b, levels = 20, color = :black, linewidth = 3.0,)

Colorbar(fig[1,3], hm, label = state_name) # , ticks = contour_levels)



for tmpax in [ax_h, ax_nh]
    tmpax.title = state_name
    tmpax.titlesize = 40
    tmpax.ylabel = "Depth [m]"
    tmpax.xlabel = "South to North [m]"
    tmpax.xlabelsize = 25
    tmpax.ylabelsize = 25 
    tmpax.limits = (extrema(y)..., extrema(z)...)
end

# ax1 = fig[jj,ii] = Axis(fig)

display(fig)

#=
ax.xticks = ([-80, -60, -30, 0, 30, 60, 80], ["80S", "60S", "30S", "0", "30N", "60N", "80N"])
pressure_levels = [1000, 850, 700, 550, 400, 250, 100, 10]
ax.yticks = (pressure_levels .* 1e2, string.(pressure_levels))
ax.yreversed = true
=#
