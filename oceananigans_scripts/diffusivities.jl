using JLD2, LinearAlgebra, Statistics

# prefix = "fluxernathy_tracers_restarted_j1_k10_averages.jld2"
# prefix = "fluxernathy_tracers_restarted_j1_k20_averages.jld2"
# prefix = "fluxernathy_tracers_restarted_smooth_forcing_j2_k2_averages.jld2"
prefix = "fluxernathy_tracers_restarted_smooth_forcing_j10_k10_averages.jld2"
include(pwd() * "/oceananigans_scripts/utils.jl")

function c_flux_gradient(jld2_file; tracer_strings=["c1", "c2", "c3", "c4"], time_index = 22)
    jl_file = jld2_file 
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
    cs   = []
    cys  = []
    czs  = []
    vcps = []
    wcps = []

    ∇c∇ᵖbs  = []  # isopycnal tracer gradients ∇c ⋅ ∇ᵖb
    ∇c∇bs  = []  # diapycnal tracer gradients ∇c ⋅ ∇b
    u⃗c∇ᵖbs = [] # isopycnal eddy flux u⃗'c' ⋅ ∇ᵖb
    u⃗c∇bs = [] # diapycnal eddy flux u⃗'c' ⋅ ∇b

    for tracer_string in tracer_strings
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
    return cs, cys, czs, vcps, wcps, ∇c∇bs, ∇c∇ᵖbs, u⃗c∇bs, u⃗c∇ᵖbs
end


## Grab Tracers for analysis
jl_file = jldopen(prefix, "r+")
cs, cys, czs, vcps, wcps, ∇c∇bs, ∇c∇ᵖbs, u⃗c∇bs, u⃗c∇ᵖbs = c_flux_gradient(jl_file)
# to incorporate multiple files, use: [cs[1:2]..., cs[3:4]...]

m,n = size(∇c∇bs[1])

K = zeros(2,2, m,n)

mask = zeros(m,n)
A = zeros(2,2)
b1 = zeros(2)
b2 = zeros(2)

numtracers = 4
case_weights = ones(numtracers) ./ numtracers # evenly weight different cases (silly due to normalization of c?)

for j in 1:m, k in 1:n
    mask[j, k] = 1.0
    for i in 1:numtracers
        ωⁱ = case_weights[i]
        # calculate matrix entries
        A[1,1] += ωⁱ * mean(∇c∇bs[i]  .* ∇c∇bs[i]  .* mask) 
        A[1,2] += ωⁱ * mean(∇c∇bs[i]  .* ∇c∇ᵖbs[i] .* mask) 
        A[2,2] += ωⁱ * mean(∇c∇ᵖbs[i] .* ∇c∇ᵖbs[i] .* mask) 
        # calculate rhs 1 for κ¹¹, κ¹²
        b1[1]  -= ωⁱ * mean(∇c∇bs[i]  .* u⃗c∇bs[i] .* mask)
        b1[2]  -= ωⁱ * mean(∇c∇ᵖbs[i] .* u⃗c∇bs[i] .* mask) 
        # calculate rhs 2 for κ²¹, κ²²
        b2[1]  -= ωⁱ * mean(∇c∇bs[i]  .* u⃗c∇ᵖbs[i] .* mask)
        b2[2]  -= ωⁱ * mean(∇c∇ᵖbs[i] .* u⃗c∇ᵖbs[i] .* mask)  
    end
    A[2,1] = A[1,2]
    # save 
    κ¹¹, κ¹² = A \ b1
    κ²¹, κ²² = A \ b2
    # save diffusivity tensor at each location in space
    K[1,1,j,k] = κ¹¹
    K[1,2,j,k] = κ¹²
    K[2,1,j,k] = κ²¹
    K[2,2,j,k] = κ²²
    # reset
    A .= 0.0
    b1 .= 0.0
    b2 .= 0.0
    mask[j,k] = 0.0
    println("finishing j = $j and k = $k")
end