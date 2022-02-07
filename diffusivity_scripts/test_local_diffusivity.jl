include(pwd() * "/oceananigans_scripts/utils.jl")
include(pwd() * "/diffusivity_scripts/utils_file.jl")

using JLD2, LinearAlgebra, Statistics

check_answer = false

case = "attempt5"
if case[1:5] == "trial"
    i = Meta.parse(case[6:end])
    jlist = [0, 1, 2]
    klist = [4 * i + 0, 4 * i + 1, 4 * i + 2, 4 * i + 3]
elseif case[1:7] == "attempt"
    i = Meta.parse(case[8:end])
    jlist = [2 * i + 3, 2 * i + 4]
    klist = [1, 2, 3, 4, 5, 6]
end

# choose file to loads
prefix = "relaxation_channel_tracers_restarted_smooth_forcing_case_" * case * "_averages.jld2"

# pick out tracer string names based on file convention
tracer_strings = []
for j in jlist, k in klist
    push!(tracer_strings, "c_j" * string(j) * "_k" * string(k))
end

## Grab Tracers for analysis
jl_file = jldopen(prefix, "r+")

cs, cys, czs, vcps, wcps, ∇c∇bs, ∇c∇ᵖbs, u⃗c∇bs, u⃗c∇ᵖbs = c_flux_gradient(jl_file, tracer_strings = tracer_strings)
# to incorporate multiple files, use: [cs[1:2]..., cs[3:4]...]
tracer_indices = 1:12
m, n = size(∇c∇bs[1])

check_K = zeros(2, 2, m, n)
y, z = get_grid(cs[1], jl_file; ghost = 3)
y̅ = avg(y)
z̅ = avg(z)
# create analytic answer, the tracer diffusivity should be the same across all tracers
for j in 1:m, k in 1:n, i in tracer_indices
    check_K[1, 1, j, k] = 3 * ((y̅[j] / y̅[end]^2)^2 + (z̅[k] / z̅[1])^2)
    check_K[1, 2, j, k] = (y̅[j] / y̅[end]^2)^2
    check_K[2, 1, j, k] = (z̅[k] / z̅[1]^2)^2
    check_K[2, 2, j, k] = 10 * ((y̅[j] / y̅[end]^2)^2 + (z̅[k] / z̅[1])^2)
    u⃗c∇bs[i][j, k] = -check_K[1, 1, j, k] * ∇c∇bs[i][j, k] - check_K[1, 2, j, k] * ∇c∇ᵖbs[i][j, k]
    u⃗c∇ᵖbs[i][j, k] = -check_K[2, 1, j, k] * ∇c∇bs[i][j, k] - check_K[2, 2, j, k] * ∇c∇ᵖbs[i][j, k]
end


K = zeros(2, 2, m, n)

mask = zeros(m, n)
A = zeros(2, 2)
b1 = zeros(2)
b2 = zeros(2)

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
        A[1, 1] += ωⁱ * mean(∇c∇bs[i] .* ∇c∇bs[i] .* mask)
        A[1, 2] += ωⁱ * mean(∇c∇bs[i] .* ∇c∇ᵖbs[i] .* mask)
        A[2, 2] += ωⁱ * mean(∇c∇ᵖbs[i] .* ∇c∇ᵖbs[i] .* mask)
        # calculate rhs 1 for κ¹¹, κ¹²
        b1[1] -= ωⁱ * mean(∇c∇bs[i] .* u⃗c∇bs[i] .* mask)
        b1[2] -= ωⁱ * mean(∇c∇ᵖbs[i] .* u⃗c∇bs[i] .* mask)
        # calculate rhs 2 for κ²¹, κ²²
        b2[1] -= ωⁱ * mean(∇c∇bs[i] .* u⃗c∇ᵖbs[i] .* mask)
        b2[2] -= ωⁱ * mean(∇c∇ᵖbs[i] .* u⃗c∇ᵖbs[i] .* mask)
    end
    A[2, 1] = A[1, 2]
    # save 
    κ¹¹, κ¹² = A \ b1
    κ²¹, κ²² = A \ b2
    # save diffusivity tensor at each location in space
    K[1, 1, j, k] = κ¹¹
    K[1, 2, j, k] = κ¹²
    K[2, 1, j, k] = κ²¹
    K[2, 2, j, k] = κ²²
    # println("The error of estimating K is ", K[:, :, j, k] - check_K[:, :, j, k])
    # reset
    A .= 0.0
    b1 .= 0.0
    b2 .= 0.0
    mask[j, k] = 0.0
    println("finishing j = $j and k = $k")
end

err = norm(check_K - K) / norm(check_K)
println("The error is $err")