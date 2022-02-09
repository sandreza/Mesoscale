using JLD2, LinearAlgebra, Statistics

include(pwd() * "/oceananigans_scripts/utils.jl")
include(pwd() * "/diffusivity_scripts/utils_file.jl")

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

tracer_indices = 1:12
m, n = size(cys[1])

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
    vcps[i][j, k] = -check_K[1, 1, j, k] * cys[i][j, k] - check_K[1, 2, j, k] * czs[i][j, k]
    wcps[i][j, k] = -check_K[2, 1, j, k] * cys[i][j, k] - check_K[2, 2, j, k] * czs[i][j, k]
end

K = zeros(2, 2, m, n)

mask = zeros(m, n)
tg = 6 # total gradients
A = zeros(tg, tg) # number of gradient locations = 1y + 1z
AA = similar(A)
B = zeros(tg, 2) # trying to predict two fluxes
b1 = zeros(2)
b2 = zeros(2)

numtracers = length(tracer_strings)

case_weights = zeros(numtracers)
for i in tracer_indices
    case_weights[i] = 1.0 / maximum(abs.(cs[i])) # normalize case by tracer magnitudes
end

j = 100
k = 10
A .= 0.0
AA .= 0.0
B .= 0.0
for i in 1:12
    ωⁱ = case_weights[i]
    # calculate matrix entries
    v⃗ = [cys[i][j, k], czs[i][j, k], czs[i][j, k+1], czs[i][j, k-1], czs[i][j-1, k], czs[i][j+1, k]]
    b = [vcps[i][j, k] wcps[i][j, k]]
    AA .+= (v⃗ * v⃗' * ωⁱ)
    B  .+= -v⃗ * [vcps[i][j, k] wcps[i][j, k]] * ωⁱ
end

# save 
KK = (AA \ B)'

i = 1
v⃗ = [cys[i][j, k], czs[i][j, k], czs[i][j, k+1], czs[i][j, k-1], czs[i][j-1, k], czs[i][j+1, k]]
b = [vcps[i][j, k] wcps[i][j, k]]
err = KK * v⃗ + b'
err[1] / abs(b[1]) * 100
err[2] / abs(b[2]) * 100
# both the following should be zero 
# the first two columns of KK check dependence on local terms 
# the subsequent columns check nonlocality
norm(check_K[:, :, j, k] - KK[1:2, 1:2])
norm(KK[1:2, 3:6])

