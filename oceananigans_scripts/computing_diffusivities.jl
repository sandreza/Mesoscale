using JLD2, LinearAlgebra, Statistics, GLMakie

# prefix = "fluxernathy_tracers_restarted_j1_k1_averages.jld2"
# prefix = "fluxernathy_nh_ns_averages.jld2"
# prefix = "fluxernathy_averages.jld2"
# prefix = "fluxernathy_tracers_restarted_j1_k10_averages.jld2"
# prefix = "fluxernathy_tracers_restarted_j1_k20_averages.jld2"
# prefix = "fluxernathy_tracers_restarted_j128_k1_averages.jld2"
# prefix = "check_h200.jld2"
# prefix = "smooth.jld2"
# prefix = "forcing_10_10.jld2"
# prefix = "fluxernathy_tracers_restarted_smooth_forcing_j2_k2_averages.jld2"
prefix = "fluxernathy_tracers_restarted_smooth_forcing_j10_k10_averages.jld2"
# prefix = "fluxernathy_tracers_restarted_smooth_forcing_j30_k3_averages.jld2"
# prefix = "fluxernathy_tracers_restarted_smooth_forcing_j6_k10_averages.jld2"
include("utils.jl")
jl_file = jldopen(prefix, "r+")



time_index = 22
b = get_field("b", time_index, jl_file) 
y,z = get_grid(b, jl_file; ghost = 3)

y̅ = avg(y)
z̅ = avg(z)
bz_p = ∂z(b, z)
by_p = ∂y(b, y)

bz = avgy(bz_p)
by = avgz(by_p)

v = get_field("v", time_index, jl_file)
w = get_field("w", time_index, jl_file)

# tracer dependent 
cs   = []
cys  = []
czs  = []
vcps = []
wcps = []

∇c∇ᵖbs  = []  # isopycnal tracer gradients ∇c ⋅ ∇ᵖb
∇c∇bs  = []  # diapycnal tracer gradients ∇c ⋅ ∇b
u⃗c∇ᵖbs = [] # isopycnal eddy flux u⃗'c' ⋅ ∇ᵖb
u⃗c∇bs = [] # diapycnal eddy flux u⃗'c' ⋅ ∇b

for tracer_string in ["c1", "c2", "c3", "c4"]
    local c = get_field(tracer_string, time_index, jl_file) 
    local vc = get_field("v"*tracer_string, time_index, jl_file)
    local wc = get_field("w"*tracer_string, time_index, jl_file)

    local vcp = avgy(avgz( avgy(vc) - avgy(v) .* c ))
    local wcp = avgy(avgz( avgz(wc) - avgz(w) .* c ))

    local cz = avgy(∂z(c, z))
    local cy = avgz(∂y(c, y))

    local vcp_avg = avgz(avgy(vcp))
    local wcp_avg = avgz(avgy(wcp))

    local ∇c∇b  = cy .* by + cz .* bz
    local ∇c∇ᵖb = cy .* bz - cz .* by

    local u⃗c∇b  = vcp .* by + wcp .* bz
    local u⃗c∇ᵖb = vcp .* bz - wcp .* by

    push!(cs, c)
    push!(cys, cy)
    push!(czs, cz)
    push!(vcps, vcp_avg)
    push!(wcps, wcp_avg)

    push!(∇c∇bs, ∇c∇b)
    push!(∇c∇ᵖbs, ∇c∇ᵖb)
    push!(u⃗c∇bs, u⃗c∇b)
    push!(u⃗c∇ᵖbs, u⃗c∇ᵖb)
end

ti = 2 # tracer index
mask = zeros(size(∇c∇bs[1]))
# mask[40:255-40, 10:30] .= 1.0

mask[40+20:255-40-20, 10+5:30-5] .= 1.0
zlevels = collect(3:31)
depth_dependent_diffusivity_estimate = Float64[]
depth_level = Float64[]
tensors = []
s_tensors = []
a_tensors = []
for zlevel in zlevels
    # mask[40+20:255-40-20, zlevel-2:zlevel+2] .= 1.0
    local ylevel = 128
    local mask = zeros(size(∇c∇bs[1]))
    local mask[ylevel-2:ylevel+2, zlevel-1:zlevel+1] .= 1.0
    # mask[40:255-40, 10:30]
    # mask[100:150, :] .= 0.0
    # mask[:,:] .= 1.0

    local ∇c∇b  = ∇c∇bs[ti]  .* mask
    local ∇c∇ᵖb = ∇c∇ᵖbs[ti] .* mask


    local u⃗c∇b  = u⃗c∇bs[ti] .* mask
    local u⃗c∇ᵖb = u⃗c∇ᵖbs[ti] .* mask

    # analytic check
    # local u⃗c∇b  = 10.0 .* ∇c∇b + -2.0 .* ∇c∇ᵖb
    # local u⃗c∇ᵖb = 1.0 .* ∇c∇b + 200.0 .* ∇c∇ᵖb


    local A = [mean(∇c∇b .* ∇c∇b) mean(∇c∇b .* ∇c∇ᵖb); mean(∇c∇b .* ∇c∇ᵖb) mean(∇c∇ᵖb .* ∇c∇ᵖb)]
    local b⃗ =  -1.0 .* [mean(∇c∇b .* u⃗c∇b) ; mean(∇c∇ᵖb .* u⃗c∇b)]
    local κ¹¹, κ¹² = A \ b⃗

    local b⃗ = -1.0 .* [mean(∇c∇b .* u⃗c∇ᵖb) ; mean(∇c∇ᵖb .* u⃗c∇ᵖb)]
    local κ²¹, κ²² = A \ b⃗

    local K = [κ¹¹ κ¹²; κ²¹ κ²²]
    local S = (K + K') * 0.5
    local A = (K - K') * 0.5

    push!(tensors, K)
    push!(s_tensors, S)
    push!(a_tensors, A)

    println("----------")
    println("zlevel is ", z[zlevel])
    local eigvalK = eigvals(K)
    println("The eigenvalues of K are ", eigvalK)
    println("---------------")
    push!(depth_dependent_diffusivity_estimate, maximum(eigvalK))
    push!(depth_level, z[zlevel])

    Ld = sqrt(mean( (u⃗c∇b .+  κ¹¹ .*  ∇c∇b .+  κ¹² .* ∇c∇ᵖb) .^2 )) / sqrt( mean( u⃗c∇b .^2 ) )
    Li = sqrt(mean( (u⃗c∇ᵖb .+  κ²¹ .*  ∇c∇b .+  κ²² .* ∇c∇ᵖb) .^2 )) / sqrt( mean(u⃗c∇ᵖb .^2 ) )
end

println("The relative standard deviation of the depth dependent diffusivity is ", sqrt(var(depth_dependent_diffusivity_estimate)) / mean(depth_dependent_diffusivity_estimate))

# For figures: 
mask = zeros(size(∇c∇bs[1]))
# mask[40:255-40, 10:30] .= 1.0

zlevel = 25
ylevel = 128
mask[ylevel-2:ylevel+2, zlevel-1:zlevel+1] .= 1.0

∇c∇b  = ∇c∇bs[ti]  .* mask
∇c∇ᵖb = ∇c∇ᵖbs[ti] .* mask

u⃗c∇b  = u⃗c∇bs[ti] .* mask
u⃗c∇ᵖb = u⃗c∇ᵖbs[ti] .* mask

A = [mean(∇c∇b .* ∇c∇b) mean(∇c∇b .* ∇c∇ᵖb); mean(∇c∇b .* ∇c∇ᵖb) mean(∇c∇ᵖb .* ∇c∇ᵖb)]
println("A for zlevel ", zlevel, " is ", A)
b⃗ =  -1.0 .* [mean(∇c∇b .* u⃗c∇b) ; mean(∇c∇ᵖb .* u⃗c∇b)]
κ¹¹, κ¹² = A \ b⃗

b⃗ = -1.0 .* [mean(∇c∇b .* u⃗c∇ᵖb) ; mean(∇c∇ᵖb .* u⃗c∇ᵖb)]
κ²¹, κ²² = A \ b⃗

K = [κ¹¹ κ¹²; κ²¹ κ²²]
S = (K + K') * 0.5
A = (K - K') * 0.5

println("----------")
println("zlevel is ", z[zlevel])
eigvalK = eigvals(K)
println("The eigenvalues of K are ", eigvalK)
println("---------------")

Ld = sqrt(mean( (u⃗c∇b .+  κ¹¹ .*  ∇c∇b .+  κ¹² .* ∇c∇ᵖb) .^2 )) / sqrt( mean( u⃗c∇b .^2 ) )
Li = sqrt(mean( (u⃗c∇ᵖb .+  κ²¹ .*  ∇c∇b .+  κ²² .* ∇c∇ᵖb) .^2 )) / sqrt( mean(u⃗c∇ᵖb .^2 ) )

lq, uq = [0.05, 0.95]
residuald = (u⃗c∇b .+  κ¹¹ .*  ∇c∇b .+  κ¹² .* ∇c∇ᵖb) ./ sqrt( mean( u⃗c∇b .^2 ) )
residuali = (u⃗c∇ᵖb .+  κ²¹ .*  ∇c∇b .+  κ²² .* ∇c∇ᵖb) ./ sqrt( mean(u⃗c∇ᵖb .^2 ) )
climsi = quantile.(Ref(residuali[:]), [lq, uq])
fig, ax, hm = heatmap(y̅, z̅, residuali, colormap = :balance, colorrange = climsi, interpolate = true)

# scatter(depth_dependent_diffusivity_estimate, zlevels)