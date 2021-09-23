using JLD2, LinearAlgebra, Statistics

# prefix = "fluxernathy_tracers_restarted_j1_k1_averages.jld2"
# prefix = "fluxernathy_nh_ns_averages.jld2"
# prefix = "fluxernathy_averages.jld2"
# prefix = "fluxernathy_tracers_restarted_j1_k10_averages.jld2"
# prefix = "fluxernathy_tracers_restarted_j1_k20_averages.jld2"
prefix = "fluxernathy_tracers_restarted_j128_k1_averages.jld2"
# prefix = "check_h200.jld2"
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

ti = 1 # tracer index
mask = zeros(size(∇c∇bs[1]))
mask[40:255-40, 10:30] .= 1.0
mask[100:150, :] .= 0.0
# mask[:,:] .= 1.0

∇c∇b  = ∇c∇bs[ti]  .* mask
∇c∇ᵖb = ∇c∇ᵖbs[ti] .* mask

u⃗c∇b  = u⃗c∇bs[ti] .* mask
u⃗c∇ᵖb = u⃗c∇ᵖbs[ti] .* mask

A = [mean(∇c∇b .* ∇c∇b) mean(∇c∇b .* ∇c∇ᵖb); mean(∇c∇b .* ∇c∇ᵖb) mean(∇c∇ᵖb .* ∇c∇ᵖb)]
b⃗ =  -1.0 .* [mean(∇c∇b .* u⃗c∇b) ; mean(∇c∇ᵖb .* u⃗c∇b)]
κ¹¹, κ¹² = A \ b⃗

b⃗ = -1.0 .* [mean(∇c∇b .* u⃗c∇ᵖb) ; mean(∇c∇ᵖb .* u⃗c∇ᵖb)]
κ²¹, κ²² = A \ b⃗

K = [κ¹¹ κ¹²; κ²¹ κ²²]
S = (K + K') * 0.5
A = (K - K') * 0.5

println("The eigenvalues of K are ", eigvals(K))

Ld = sqrt(mean( (u⃗c∇b .+  κ¹¹ .*  ∇c∇b .+  κ¹² .* ∇c∇ᵖb) .^2 )) / sqrt( mean( u⃗c∇b .^2 ) )
Li = sqrt(mean( (u⃗c∇ᵖb .+  κ²¹ .*  ∇c∇b .+  κ²² .* ∇c∇ᵖb) .^2 )) / sqrt( mean(u⃗c∇ᵖb .^2 ) )

lq, uq = [0.05, 0.95]
residuald = (u⃗c∇b .+  κ¹¹ .*  ∇c∇b .+  κ¹² .* ∇c∇ᵖb) ./ sqrt( mean( u⃗c∇b .^2 ) )
residuali = (u⃗c∇ᵖb .+  κ²¹ .*  ∇c∇b .+  κ²² .* ∇c∇ᵖb) ./ sqrt( mean(u⃗c∇ᵖb .^2 ) )
climsi = quantile.(Ref(residuali[:]), [lq, uq])
fig, ax, hm = heatmap(y̅, z̅, residuali, colormap = :balance, colorrange = climsi, interpolate = true)