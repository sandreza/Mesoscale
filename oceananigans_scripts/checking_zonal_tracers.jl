using JLD2
using GLMakie
using Statistics

tracer_case_j = 1
tracer_case_k = 1
# file output
# prefix = "fluxernathy_tracers_iteration9986400.jld2"
# prefix = "fluxernathy_tracers_restarted_j1_k1_averages.jld2"
# prefix = "checkme.jld2"
# prefix = "fluxernathy_nh_ns_averages.jld2"
# prefix = "fluxernathy_averages.jld2"
prefix = "fluxernathy_tracers_restarted_j1_k10_averages.jld2"
jl_file = jldopen(prefix, "r+")

function get_grid(field, jl_file; ghost = 3)

    yC = collect(jl_file["grid"]["yC"])[ghost+1:end-ghost]
    zC = collect(jl_file["grid"]["zC"])[ghost+1:end-ghost]

    yF = collect(jl_file["grid"]["yF"])[ghost+1:end-ghost]
    zF = collect(jl_file["grid"]["zF"])[ghost+1:end-ghost]


    y = size(field)[1] == size(yC)[1] ? yC : yF
    z = size(field)[2] == size(zC)[1] ? zC : zF

    return y, z
end

function get_grid2(field, jl_file; ghost = 3)

    yC = collect(jl_file["grid"]["yᵃᶜᵃ"])[ghost+1:end-ghost]
    zC = collect(jl_file["grid"]["zᵃᵃᶜ"])[ghost+1:end-ghost]

    yF = collect(jl_file["grid"]["yᵃᶠᵃ"])[ghost+1:end-ghost]
    zF = collect(jl_file["grid"]["zᵃᵃᶠ"])[ghost+1:end-ghost]


    y = size(field)[1] == size(yC)[1] ? yC : yF
    z = size(field)[2] == size(zC)[1] ? zC : zF

    return y, z
end

function get_field(field_string::String, time_index::Int, jl_file::JLD2.JLDFile; ghost = 3)
    t_keys = keys(jl_file["timeseries"][field_string])
    field = jl_file["timeseries"][field_string][t_keys[time_index]]
    return field[1,:,:]
end

function ∂y(field, y)
    δf = field[2:end, :] - field[1:end-1, :]
    δy = y[2:end] - y[1:end-1]
    δy = reshape(δy[:], (length(δy), 1))
    return δf ./ δy
end

function ∂z(field,z)
    δf = field[:, 2:end] - field[:, 1:end-1]
    δz = z[2:end] - z[1:end-1]
    δz = reshape(δz[:], (1, length(δz)))
    return δf ./ δz
end

function avg(field)
    return (field[2:end] + field[1:end-1]) .* 0.5
end
function avgy(field)
    return (field[2:end, :] + field[1:end-1, :]) .* 0.5
end
function avgz(field)
    return (field[:, 2:end] + field[:, 1:end-1]) .* 0.5
end

#=
field = get_field("c1", 9, jl_file)
y, z = get_grid(field, jl_file)

y̅ = avg(y)
z̅ = avg(z)
fieldz = ∂z(field,z) 
# heatmap(y, z, field, colormap = :balance, interpolate = true)

heatmap(y, z̅, fieldz, colormap = :balance, interpolate = true)


field1 = get_field("c3", 9, jl_file)
field2 = get_field("c4", 9, jl_file)
field1_m = ∂z(field1, z)
field2_m = ∂z(field2, z)

fig, ax, hm = heatmap(y, z̅, field1_m, colormap = :balance, interpolate = true)
ax2 = fig[1,2] = Axis(fig)
heatmap!(ax2, y, z̅, field2_m, colormap = :balance, interpolate = true)

ax3 = fig[1,3] = Axis(fig)
cc = maximum(abs.(field1_m))
clims = (-cc, cc)
heatmap!(ax3, y, z̅, field1_m - field2_m, colormap = :balance, interpolate = true, colorrange = clims)
=#

# data pairs
#=
custom_clims = false
tracer_string = "c1"
time_index = 20
b = get_field("b", time_index, jl_file) 
c = get_field(tracer_string, time_index, jl_file) 
v = get_field("v", time_index, jl_file)
w = get_field("w", time_index, jl_file)
vc = get_field("v"*tracer_string, time_index, jl_file)
wc = get_field("w"*tracer_string, time_index, jl_file)
vcp = avgy(vc) - avgy(v) .* c
wcp = avgz(wc) - avgz(w) .* c

uq = 0.95 #upper quantile
lq = 0.05 #lower quantile

y,z = get_grid(c, jl_file; ghost = 3)

y̅ = avg(y)
z̅ = avg(z)

cz = ∂z(c, z)
cy = ∂y(c, y)

bz = ∂z(b, z)
by = ∂y(b, y)

if custom_clims
    clims = (0, 2e-5)
else
    clims = quantile.(Ref(cz[:]), [lq, uq])
end

fig, ax, hm = heatmap(y, z̅, cz, colormap = :balance, interpolate = true, colorrange = clims)
ax.title = "∂z(" * tracer_string * ")"
Colorbar(fig[1,2], hm)

if custom_clims
    clims2 = (-2e-9, 1e-8)
else
    clims2 = quantile.(Ref(cy[:]), [lq, uq])
end

ax2 = fig[1,3] = Axis(fig)
hm2 = heatmap!(ax2, y̅, z, cy, colormap = :balance, interpolate = true, colorrange = clims2)
ax2.title = "∂y(" * tracer_string * ")"
Colorbar(fig[1,4], hm2)

if custom_clims
    clims3 = (-1e-9, 2e-8)
else
    clims3 = quantile.(Ref(wcp[:]), [lq,uq])
end

ax3 = fig[2,1] = Axis(fig)
hm3 = heatmap!(ax3, y, z, wcp, colormap = :balance, interpolate = true, colorrange = clims3)
ax3.title = "w'"* tracer_string * "'"
Colorbar(fig[2,2], hm3)

if custom_clims
    clims4 = (-2e-5, 0)
else
    clims4 = quantile.(Ref(vcp[:]), [lq,uq])
end

ax4 = fig[2,3] = Axis(fig)
hm4 = heatmap!(ax4, y, z, vcp, colormap = :balance, interpolate = true, colorrange = clims4)
ax4.title = "v'"* tracer_string * "'"
Colorbar(fig[2,4], hm4)

##
vcp_avg = avgz(avgy(vcp))
wcp_avg = avgz(avgy(wcp))
bz_avg = avgy(bz)
by_avg = avgz(by)

dcp = vcp_avg .* by_avg + wcp_avg .* bz_avg
icp = vcp_avg .* bz_avg - wcp_avg .* by_avg

if custom_clims
    clims5 = quantile.(Ref(dcp[:]), [lq,uq])
else
    clims5 = quantile.(Ref(dcp[:]), [lq,uq])
end

if custom_clims
    clims6 = quantile.(Ref(icp[:]), [lq,uq])
    clims5 = clims6
else
    clims6 = quantile.(Ref(icp[:]), [lq,uq])
    clims5 = clims6
end

ax5 = fig[3,1] = Axis(fig)
hm5 = heatmap!(ax5, y̅, z̅, dcp , colormap = :balance, interpolate = true, colorrange = clims5)
ax5.title = "u⃗'"* tracer_string * "'" * "⋅∇b"
Colorbar(fig[3,2], hm5)

ax6 = fig[3,3] = Axis(fig)
hm6 = heatmap!(ax6, y̅, z̅, icp, colormap = :balance, interpolate = true, colorrange = clims6)
ax6.title = "u⃗'"* tracer_string * "'" * "⋅∇ᵖb"
Colorbar(fig[3,4], hm6)
=#


clim_options = [:custom, :symmetric]
custom_clims =  clim_options[2] 
tracer_string = "c4"
time_index = 20
b = get_field("b", time_index, jl_file) 
c = get_field(tracer_string, time_index, jl_file) # -get_field("c3", time_index, jl_file) 
v = get_field("v", time_index, jl_file)
w = get_field("w", time_index, jl_file)
vc = get_field("v"*tracer_string, time_index, jl_file)
wc = get_field("w"*tracer_string, time_index, jl_file)
vcp = avgy(vc) - avgy(v) .* c
wcp = avgz(wc) - avgz(w) .* c

uq = 0.95 #upper quantile
lq = 0.05 #lower quantile

y,z = get_grid(c, jl_file; ghost = 3)

y̅ = avg(y)
z̅ = avg(z)

cz = ∂z(c, z)
cy = ∂y(c, y)

bz = ∂z(b, z)
by = ∂y(b, y)

vcp_avg = avgz(avgy(vcp))
wcp_avg = avgz(avgy(wcp))
bz_avg = avgy(bz)
by_avg = avgz(by)



fig = Figure()


clims01 = quantile.(Ref(c[:]), [lq, uq])
ax01 = fig[1,1]= Axis(fig)
hm01 = heatmap!(ax01, y, z, c, colormap = :thermometer, interpolate = true, colorrange = clims01)
ax01.title = tracer_string 
Colorbar(fig[1,2], hm01)


∇b = sqrt.(bz_avg .^2 + by_avg .^2)
clims03 = quantile.(Ref(∇b[:]), [lq, uq])
ax03 = fig[3,1]= Axis(fig)
hm03 = heatmap!(ax03, y, z, ∇b, colormap = :thermometer, interpolate = true, colorrange = clims03)
ax03.title = "|∇b|"
Colorbar(fig[3,2], hm03)



if custom_clims==:custom  
    clims = (0, 2e-5)
elseif custom_clims==:symmetric
    clims = quantile.(Ref(cz[:]), [lq, uq])
    cu = maximum(abs.(clims))
    clims = (-cu, cu)
else
    clims = quantile.(Ref(cz[:]), [lq, uq])
end

ax = fig[1,1+2]= Axis(fig)
hm = heatmap!(ax, y, z̅, cz, colormap = :balance, interpolate = true, colorrange = clims)
ax.title = "∂z(" * tracer_string * ")"
Colorbar(fig[1,2+2], hm)

if custom_clims==:custom 
    clims2 = (-2e-9, 1e-8)
elseif custom_clims==:symmetric
    clims2 = quantile.(Ref(cy[:]), [lq, uq])
    cu = maximum(abs.(clims2))
    clims2 = (-cu, cu)
else
    clims2 = quantile.(Ref(cy[:]), [lq, uq])
end

ax2 = fig[1,3+2] = Axis(fig)
hm2 = heatmap!(ax2, y̅, z, cy, colormap = :balance, interpolate = true, colorrange = clims2)
ax2.title = "∂y(" * tracer_string * ")"
Colorbar(fig[1,4+2], hm2)

if custom_clims==:custom
    clims3 = (-1e-9, 2e-8)
elseif custom_clims==:symmetric
    clims3 = quantile.(Ref(wcp[:]), [lq, uq])
    cu = maximum(abs.(clims3))
    clims3 = (-cu, cu)
else
    clims3 = quantile.(Ref(wcp[:]), [lq,uq])
end

ax3 = fig[2,1+2] = Axis(fig)
hm3 = heatmap!(ax3, y, z, wcp, colormap = :balance, interpolate = true, colorrange = clims3)
ax3.title = "w'"* tracer_string * "'"
Colorbar(fig[2,2+2], hm3)

if custom_clims==:custom
    clims4 = (-2e-5, 0)
elseif custom_clims==:symmetric
    clims4 = quantile.(Ref(vcp[:]), [lq, uq])
    cu = maximum(abs.(clims4))
    clims4 = (-cu, cu)
else
    clims4 = quantile.(Ref(vcp[:]), [lq,uq])
end

ax4 = fig[2,3+2] = Axis(fig)
hm4 = heatmap!(ax4, y, z, vcp, colormap = :balance, interpolate = true, colorrange = clims4)
ax4.title = "v'"* tracer_string * "'"
Colorbar(fig[2,4+2], hm4)

##
vcp_avg = avgz(avgy(vcp))
wcp_avg = avgz(avgy(wcp))
bz_avg = avgy(bz)
by_avg = avgz(by)

dcp = vcp_avg .* by_avg + wcp_avg .* bz_avg
icp = vcp_avg .* bz_avg - wcp_avg .* by_avg

if custom_clims==:custom
    clims5 = quantile.(Ref(dcp[:]), [lq,uq])
elseif custom_clims==:symmetric
    clims5 = quantile.(Ref(dcp[:]), [lq, uq])
    cu = maximum(abs.(clims5))
    clims5 = (-cu, cu)
else
    clims5 = quantile.(Ref(dcp[:]), [lq,uq])
end

if custom_clims==:custom
    clims6 = quantile.(Ref(icp[:]), [lq,uq])
    clims5 = clims6
elseif custom_clims==:symmetric
    clims6 = quantile.(Ref(icp[:]), [lq, uq])
    cu = maximum(abs.(clims6))
    clims6 = (-cu, cu)
    # clims5 = clims6
else
    clims6 = quantile.(Ref(icp[:]), [lq,uq])
    clims5 = clims6
end

ax5 = fig[3,1+2] = Axis(fig)
hm5 = heatmap!(ax5, y̅, z̅, dcp , colormap = :balance, interpolate = true, colorrange = clims5)
ax5.title = "u⃗'"* tracer_string * "'" * "⋅∇b"
Colorbar(fig[3,2+2], hm5)

ax6 = fig[3,3+2] = Axis(fig)
hm6 = heatmap!(ax6, y̅, z̅, icp, colormap = :balance, interpolate = true, colorrange = clims6)
ax6.title = "u⃗'"* tracer_string * "'" * "⋅∇ᵖb"
Colorbar(fig[3,4+2], hm6)

