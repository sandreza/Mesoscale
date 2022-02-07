using JLD2, LinearAlgebra, Statistics
plotting = false
check_answer = false
# prefix = "fluxernathy_tracers_restarted_j1_k10_averages.jld2"
# prefix = "fluxernathy_tracers_restarted_j1_k20_averages.jld2"
# prefix = "fluxernathy_tracers_restarted_smooth_forcing_j2_k2_averages.jld2"
# prefix = "fluxernathy_tracers_restarted_smooth_forcing_j10_k10_averages.jld2"
# prefix = "fluxernathy_tracers_restarted_smooth_forcing_j6_k10_averages.jld2"
# prefix = "fluxernathy_tracers_restarted_smooth_forcing_j6_k6_averages.jld2"
# prefix = "relaxation_channel_tracers_restarted_smooth_forcing_j1_k1_averages.jld2"
# prefix = "relaxation_channel_tracers_restarted_smooth_forcing_j1_k10_averages.jld2"
# prefix = "relaxation_channel_tracers_restarted_smooth_forcing_j10_k1_averages.jld2"
prefix = "relaxation_channel_tracers_restarted_smooth_forcing_j30_k1_averages.jld2"
prefix = "relaxation_channel_tracers_restarted_smooth_forcing_case_trial0_averages.jld2"

# case = "trial5"
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

tracer_strings = []
for j in jlist, k in klist
    push!(tracer_strings, "c_j" * string(j) * "_k" * string(k))
end
# selection of tracers later

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

#=
case = "trial0"
if case[1:5] == "trial"
    i = Meta.parse(case[6:end])
    jlist = [0,1,2]
    klist = [4*i+0,4*i+1,4*i+2,4*i+3]
end
tracer_strings = []
for j in jlist, k in klist
    push!(tracer_strings, "c_j"*string(j) * "_k"*string(k))
end
=#

cs, cys, czs, vcps, wcps, ∇c∇bs, ∇c∇ᵖbs, u⃗c∇bs, u⃗c∇ᵖbs = c_flux_gradient(jl_file, tracer_strings = tracer_strings)
# to incorporate multiple files, use: [cs[1:2]..., cs[3:4]...]

m, n = size(∇c∇bs[1])
if check_answer
    check_K = zeros(2,2, m,n)
    y,z = get_grid(cs[1], jl_file; ghost = 3)
    y̅ = avg(y)
    z̅ = avg(z)
    for j in 1:m, k in 1:n, i in 1:4
        check_K[1,1,j,k] = 3*((y̅[j]/y̅[end]^2)^2 + (z̅[k]/z̅[1])^2)
        check_K[1,2,j,k] = (y̅[j]/y̅[end]^2)^2
        check_K[2,1,j,k] = (z̅[k]/z̅[1]^2)^2
        check_K[2,2,j,k] = 10*((y̅[j]/y̅[end]^2)^2 + (z̅[k]/z̅[1])^2)
        u⃗c∇bs[i][j,k]  = -check_K[1,1,j,k]*∇c∇bs[i][j,k] - check_K[1,2,j,k] * ∇c∇ᵖbs[i][j,k]
        u⃗c∇ᵖbs[i][j,k] = - check_K[2,1,j,k]* ∇c∇bs[i][j,k] - check_K[2,2,j,k] * ∇c∇ᵖbs[i][j,k]
    end
end

K = zeros(2,2, m,n)

mask = zeros(m,n)
A = zeros(2,2)
b1 = zeros(2)
b2 = zeros(2)

tracer_indices = 1:12
numtracers = length(tracer_strings)

case_weights = zeros(numtracers)
for i in tracer_indices
    case_weights[i] = 1.0 / maximum(abs.(cs[i])) # normalize case by tracer magnitudes
end

for j in 1:m, k in 1:n
    mask[j, k] = 1.0
    for i in tracer_indices
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

if check_answer 
    err = norm(check_K - K) / norm(check_K)
    println("The error is $err")
end

if plotting
    using GLMakie 
    interpolation_flag = true
    fig = Figure()
    ax11 = fig[1,1] = Axis(fig, title = "K¹¹ : diapycnal gradient, diapycnal flux")
    field = K[1,1,:,:];
    clims = quantile.(Ref(field[:]), [0.1, 0.9])
    hm11 = heatmap!(ax11, field, colorrange = clims, colormap = :thermal, interpolate = interpolation_flag)

    Colorbar(fig[1,2], hm11, label = " ",
    topspinevisible = true, 
    bottomspinevisible = true, 
    leftspinevisible = true,
    rightspinevisible = true)

    ax12 = fig[1,3] = Axis(fig, title = "K¹² : isopycnal gradient, diapycnal flux")
    field = K[1,2,:,:];
    clims = quantile.(Ref(field[:]), [0.1, 0.9])
    hm12 = heatmap!(ax12, field, colorrange = clims, colormap = :thermal, interpolate = interpolation_flag)

    Colorbar(fig[1,4], hm12, label = " ",
    topspinevisible = true, 
    bottomspinevisible = true, 
    leftspinevisible = true,
    rightspinevisible = true)

    ax21 = fig[2,1] = Axis(fig, title = "K²¹: diapycnal gradient, isopycnal flux")
    field = K[2,1,:,:];
    clims = quantile.(Ref(field[:]), [0.1, 0.9])
    hm21 = heatmap!(ax21, field, colorrange = clims, colormap = :thermal, interpolate = interpolation_flag)

    Colorbar(fig[2,2], hm21, label = " ",
    topspinevisible = true, 
    bottomspinevisible = true, 
    leftspinevisible = true,
    rightspinevisible = true)

    ax22 = fig[2,3] = Axis(fig, title = "K²² : isopycnal gradient, isopycnal flux")
    field = K[2,2,:,:];
    clims = quantile.(Ref(field[:]), [0.1, 0.9])
    hm22 = heatmap!(ax22, field, colorrange = clims, colormap = :thermal, interpolate = interpolation_flag)

    Colorbar(fig[2,4], hm22, label = " ",
    topspinevisible = true, 
    bottomspinevisible = true, 
    leftspinevisible = true,
    rightspinevisible = true)

end
