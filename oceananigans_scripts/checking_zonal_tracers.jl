using JLD2
using GLMakie
using Statistics

case = "trial5"
case = "attempt5"
if case[1:5] == "trial"
    i = Meta.parse(case[6:end])
    jlist = [0, 1, 2]
    klist = [4*i+0,4*i+1,4*i+2,4*i+3]
elseif case[1:7] == "attempt"
    i = Meta.parse(case[8:end])
    jlist = [2 * i + 3, 2*i + 4]
    klist = [1, 2, 3, 4, 5, 6]
end

prefix = "relaxation_channel_tracers_restarted_smooth_forcing_case_"*case * "_averages.jld2"
tind = [2, 6]
tracer_string = "b"
tracer_string = "c_j"*string(jlist[tind[1]]) * "_k"*string(klist[tind[2]])
tracer_case_j = 1
tracer_case_k = 1
# file output
# prefix = "fluxernathy_tracers_iteration9986400.jld2"
# prefix = "fluxernathy_tracers_restarted_j1_k1_averages.jld2"
# prefix = "checkme2.jld2"
# prefix = "fluxernathy_nh_ns_averages.jld2"
# prefix = "fluxernathy_averages.jld2"
# prefix = "fluxernathy_tracers_restarted_j1_k10_averages.jld2"
# prefix = "fluxernathy_tracers_restarted_j1_k20_averages.jld2"
# prefix = "check_hhh.jld2"
# prefix = "fluxernathy_tracers_restarted_j128_k1_averages.jld2"
# prefix = "check_h200.jld2"
# prefix = "smooth.jld2"
# prefix = "momodes.jld2"
# prefix = "plzwork.jdl2"
# prefix = "forcing_10_10.jld2"
# prefix = "fluxernathy_tracers_restarted_smooth_forcing_j10_k10_averages.jld2"
# prefix = "fluxernathy_tracers_restarted_smooth_forcing_j30_k3_averages.jld2"
# prefix = "fluxernathy_tracers_restarted_smooth_forcing_j6_k10_averages.jld2"
# prefix = "fluxernathy_tracers_restarted_smooth_forcing_j6_k6_averages.jld2"
# prefix = "relaxation_channel_nh_averages.jld2"
# prefix = "relaxation_channel_tracers_restarted_smooth_forcing_j1_k1_averages.jld2"
# prefix = "relaxation_channel_tracers_restarted_smooth_forcing_j2_k2_averages.jld2"
# prefix = "relaxation_channel_tracers_restarted_smooth_forcing_j1_k10_averages.jld2"
# prefix = "relaxation_channel_tracers_restarted_smooth_forcing_j10_k1_averages.jld2"
# prefix = "relaxation_channel_tracers_restarted_smooth_forcing_j30_k1_averages.jld2"
# prefix = "relaxation_channel_tracers_restarted_smooth_forcing_j1_k32_averages.jld2"
# prefix = "sin_relaxation_channel_nh_averages.jld2"

# prefix = "relaxation_channel_tracers_restarted_smooth_forcing_case_trial0_averages.jld2"
jl_file = jldopen(prefix, "r+")

function get_grid(field, jl_file; ghost = 3)
    gridstrings = keys(jl_file["grid"])

    yC_string = "yC" in gridstrings ? "yC" : "yᵃᶜᵃ" 
    zC_string = "zC" in gridstrings ? "zC" : "zᵃᵃᶜ"
    yF_string = "yF" in gridstrings ? "yF" : "Δyᵃᶠᵃ"
    zF_string = "zF" in gridstrings ? "zF" : "zᵃᵃᶠ"
    yC = collect(jl_file["grid"][yC_string])[ghost+1:end-ghost]
    zC = collect(jl_file["grid"][zC_string])[ghost+1:end-ghost]

    yF = collect(jl_file["grid"][yF_string])[ghost+1:end-ghost]
    zF = collect(jl_file["grid"][zF_string])[ghost+1:end-ghost]


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



# start here
clim_options = [:custom, :symmetric]
eke_avail = true # change this nob
diff_field = false
buoyancy_contour = true
check_grad = true # change this nob
custom_clims =  clim_options[2] 
tracer_string2 = nothing # "c2" # "c3" #nothing # "c2" # nothing # "c2" # nothing
time_index = 22
b = get_field("b", time_index, jl_file) 
if tracer_string2==nothing
    c = get_field(tracer_string, time_index, jl_file) 
else
    c = get_field(tracer_string, time_index, jl_file) - get_field(tracer_string2, time_index, jl_file)
end
v = get_field("v", time_index, jl_file)
w = get_field("w", time_index, jl_file)
if tracer_string2 == nothing
    vc = get_field("v"*tracer_string, time_index, jl_file)
    wc = get_field("w"*tracer_string, time_index, jl_file)
else 
    vc = get_field("v" * tracer_string, time_index, jl_file) - get_field("v"*tracer_string2, time_index, jl_file)
    wc = get_field("w" * tracer_string, time_index, jl_file) - get_field("w"*tracer_string2, time_index, jl_file)
end


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

fig = Figure(resolution = (2800, 1400))

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

if eke_avail
    u = get_field("u", time_index, jl_file)
    v = get_field("v", time_index, jl_file)
    w = get_field("w", time_index, jl_file)
    uu = get_field("uu", time_index, jl_file)
    vv = get_field("vv", time_index, jl_file)
    ww = get_field("ww", time_index, jl_file)

    eke = uu .- (u .* u) + avgy(vv .- (v .* v) )  + avgz(ww .- (w .* w) )

    clims02 = quantile.(Ref(eke[:]), [lq, uq])
    ax02 = fig[2,1]= Axis(fig)
    hm02 = heatmap!(ax02, y, z, eke, colormap = :thermometer, interpolate = true, colorrange = clims02)
    ax02.title = "eke"
    Colorbar(fig[2,2], hm03)
else
    nothing
end

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

if buoyancy_contour
    blevels = b[end, 1:2:end]
    for axval in [ax01, ax02, ax03, ax, ax2, ax3, ax4, ax5, ax6]
        contour!(axval, y,z, b, levels = blevels, color = :black, linewidth = 1)
    end
end

if check_grad 

    ∇c∇b  = avgz(cy) .* by_avg  + avgy(cz) .* bz_avg 
    ∇c∇ᵖb = avgz(cy) .* bz_avg  - avgy(cz) .* by_avg
    
    if custom_clims==:custom
        clims7 = quantile.(Ref(∇c∇b[:]), [lq,uq])
    elseif custom_clims==:symmetric
        clims7 = quantile.(Ref(∇c∇b[:]), [lq, uq])
        cu = maximum(abs.(clims7))
        clims7 = (-cu, cu)
    else
        clims7 = quantile.(Ref(∇c∇b[:]), [lq,uq])
    end

    if custom_clims==:custom
        clims8 = quantile.(Ref(∇c∇ᵖb[:]), [lq,uq])
        clims7 = clims8
    elseif custom_clims==:symmetric
        clims8 = quantile.(Ref(∇c∇ᵖb[:]), [lq, uq])
        cu = maximum(abs.(clims8))
        clims8 = (-cu, cu)
        # clims5 = clims6
    else
        clims8 = quantile.(Ref(∇c∇ᵖb[:]), [lq,uq])
        clims7 = clims8
    end

    ax7 = fig[4,1+2+2] = Axis(fig)
    hm7 = heatmap!(ax7, y̅, z̅, -1 .* ∇c∇b , colormap = :balance, interpolate = true, colorrange = clims7)
    ax7.title = "-∇"* tracer_string * "⋅∇b = " * "-∇ᵖ"* tracer_string * "⋅∇ᵖb"
    Colorbar(fig[4,2+2+2], hm7)

    ax8 = fig[4,3+2-2] = Axis(fig)
    hm8 = heatmap!(ax8, y̅, z̅, ∇c∇ᵖb, colormap = :balance, interpolate = true, colorrange = clims8)
    ax8.title = "∇"* tracer_string * "⋅∇ᵖb = " * "-∇ᵖ"* tracer_string * "⋅∇b"
    Colorbar(fig[4,4+2-2], hm8)

    blevels = b[end, 1:2:end]

    # quick check 

    ∇c = sqrt.(avgz(cy) .* avgz(cy)  + avgy(cz) .* avgy(cz))
    clims04 = quantile.(Ref( ∇c[:]), [lq, uq])
    ax04 = fig[4,1]= Axis(fig)
    hm04 = heatmap!(ax04, y̅, z̅, ∇c , colormap = :thermometer, interpolate = true, colorrange = clims04)
    ax04.title = "|∇"* tracer_string * "|"
    Colorbar(fig[4,2], hm04)

    for axval in [ax7, ax8, ax04]
        contour!(axval, y,z, b, levels = blevels, color = :black, linewidth = 1)
    end

    for axval in [ax7, ax8, ax04]
        axval.titlesize = 40
    end
end

for axval in [ax01, ax02, ax03, ax, ax2, ax3, ax4, ax5, ax6]
    axval.titlesize = 40
end

#=
fig, ax, hm01 = heatmap(y, z, c - oldc, colormap = :thermometer, interpolate = true, colorrange = clims01)
ax.title = tracer_string 
Colorbar(fig[1,2], hm01)
=#
#=
clims6 = quantile.(Ref(w[:]), [lq,uq])
clims6 = (-clims6[2], clims6[2])
ax6 = fig[4,1] = Axis(fig)
hm6 = heatmap!(ax6, y, z, avgz(w) , colormap = :balance, interpolate = true, colorrange = clims6)
ax6.title = "w"
Colorbar(fig[4,2], hm6)


clims7 = quantile.(Ref(u[:]), [lq,uq])
ax7 = fig[4,1+2+2] = Axis(fig)
hm7 = heatmap!(ax7, y, z, u , colormap = :thermometer, interpolate = true, colorrange = clims7)
ax7.title = "u"
Colorbar(fig[4,2+2+2], hm7)

clims8 = (-1e-4, 1e-4)
ax8 = fig[4,3+2-2] = Axis(fig)
hm8 = heatmap!(ax8, y, z, avgy(v), colormap = :balance, interpolate = true, colorrange = clims8)
ax8.title = "v"
Colorbar(fig[4,4+2-2], hm8)

blevels = b[end, 1:2:end]
=#
# quick check 
#=
clims04 = quantile.(Ref( w[:]), [lq, uq])
clims04 = quantile.(Ref(w[:]), [lq, uq])
cu = maximum(abs.(clims04))
clims04 = (-cu, cu)
ax04 = fig[4,1]= Axis(fig)
hm04 = heatmap!(ax04, y, z, avgz(w) , colormap = :balance, interpolate = true, colorrange = clims04)
ax04.title = "w"
Colorbar(fig[4,2], hm04)

for axval in [ax7, ax8, ax04]
    contour!(axval, y,z, b, levels = blevels, color = :black, linewidth = 1)
end

for axval in [ax7, ax8, ax04]
    axval.titlesize = 40
end
=#

#=
# checking convergence
meanbs = Float64[]
for i in 3:30
    bb = get_field("b", i, jl_file) 
    push!(meanbs, mean(bb))
end
scatter(meanbs)

=#


#=
Ly = 2e6
windstress = @. exp( -(y- Ly/2)^2 / ( Ly^2 / 32) ) - exp( -(0- Ly/2)^2 / ( Ly^2 / 32) )
windstress2 = @. sin(π * y / Ly)
fig, ax, sc = scatter(y,  windstress)
scatter!(y, windstress2)
=#



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