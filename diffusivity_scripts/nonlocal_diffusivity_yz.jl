using JLD2, LinearAlgebra, Statistics

include(pwd() * "/oceananigans_scripts/utils.jl")
include(pwd() * "/diffusivity_scripts/utils_file.jl")
include(pwd() * "/diffusivity_scripts/spectral_bisection_utils.jl")

error_vc = []
error_wc = []

# trial0 and attempt8
diffusivities = []
cases = ["trial1", "trial2", "trial3", "trial4"]

cases = []

for i in 1:5
    push!(cases, "trial" * string(i))
end

for i in 1:7
    push!(cases, "attempt" * string(i))
end

for i in 9:16
    push!(cases, "attempt" * string(i))
end

cs = []
cys = []
czs = []
vcps = []
wcps = []
tracer_string_list = []

for case in cases
    println("currently on case ", case)
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

    l_cs, l_cys, l_czs, l_vcps, l_wcps, l_∇c∇bs, l_∇c∇ᵖbs, l_u⃗c∇bs, l_u⃗c∇ᵖbs = c_flux_gradient(jl_file, tracer_strings = tracer_strings)
    # to incorporate multiple files, use: [cs[1:2]..., cs[3:4]...]
    for i in eachindex(l_cys)
        push!(cs, l_cs[i])
        push!(cys, l_cys[i])
        push!(czs, l_czs[i])
        push!(vcps, l_vcps[i])
        push!(wcps, l_wcps[i])
        push!(tracer_string_list, tracer_strings)
    end
end

##
skip = 12 # number of tracers used to estimate diffusivity
skipper = floor(Int, length(cs) / skip)
tracer_index_partitions = [skip*(j-1)+1:skip*j for j in 1:skipper]

for tracer_indices in tracer_index_partitions
    println("currently on ", tracer_indices)
    m, n = size(cys[1])

    jn_u = 0 #jn_u, j "neighbors" right
    jn_d = 0 #jn_d, j "neighbors" left
    kn_u = 0 #kn_u, k "neighbors" up
    kn_d = 0 #kn_d, k "neighbors" down

    numtracers = length(cs)

    case_weights = zeros(numtracers)
    for i in tracer_indices
        case_weights[i] = 1.0 / maximum(abs.(cs[i])) # normalize case by tracer magnitudes
    end

    tg = 2 * (jn_u + jn_d + 1) * (kn_u + kn_d + 1)
    K = zeros(2, tg, m, n)


    @assert tg < (tracer_indices[end] - tracer_indices[1])

    for j in 1:m, k in 1:n
        # fix for boundaries
        l_jn_d = (j - jn_d < 1) ? j - 1 : jn_d
        l_kn_d = (k - kn_d < 1) ? k - 1 : kn_d
        l_jn_u = (j + jn_u > m) ? m - j : jn_u
        l_kn_u = (k + kn_u > n) ? n - k : kn_u

        ntg = 2 * (l_jn_u + l_jn_d + 1) * (l_kn_u + l_kn_d + 1)
        AA = zeros(ntg, ntg)
        BB = zeros(ntg, 2)


        for i in tracer_indices
            ωⁱ = case_weights[i]
            # calculate matrix entries,  czs[i][j, k+1], czs[i][j, k-1], czs[i][j-1, k], czs[i][j+1, k]
            v⃗ = [cys[i][j-l_jn_d:j+l_jn_u, k-l_kn_d:k+l_kn_u][:]; czs[i][j-l_jn_d:j+l_jn_u, k-l_kn_d:k+l_kn_u][:]]
            b = [vcps[i][j, k] wcps[i][j, k]]
            AA .+= (v⃗ * v⃗' * ωⁱ)
            BB .+= -v⃗ * [vcps[i][j, k] wcps[i][j, k]] * ωⁱ
        end
        # regularize coefficients
        Δ = similar(AA) .* 0
        for i in 1:ntg-1
            Δ[i, i] = 2
            Δ[i, i+1] = -1
            Δ[i+1, i] = -1
        end
        #=
        hntg = (l_jn_u + l_jn_d + 1) * (l_kn_u + l_kn_d + 1)
        # Handle boundaries (only valid for 1D)
        Δ[1, 1] = 1
        Δ[hntg, hntg] = 1
        Δ[hntg, hntg+1] = 0
        Δ[hntg, hntg-1] = -1

        Δ[hntg+1, hntg+1] = 1
        Δ[hntg+1, hntg+1+1] = -1
        Δ[hntg+1, hntg-1+1] = 0
        Δ[end, end] = 1

        λ = 1e-8 * maximum(abs.(AA))
        AA += λ * Δ
        =#
        # save
        K[:, 1:ntg, j, k] .= Transpose(AA \ BB)
        if (j == 100) & (k == 10)
            println("------------")
            println("j = ", j, " k = ", k)
            println(K[:, 1:ntg, j, k])
            # println(Transpose(AA \ BB))
            println("------------")
        end

    end
    push!(diffusivities, K)

    # Compute Error on a Per Tracer Basis
    for i in tracer_indices
        # v⃗ = [cys[i][j, k], czs[i][j, k], czs[i][j, k+1], czs[i][j, k-1], czs[i][j-1, k], czs[i][j+1, k]]
        # b = [vcps[i][j, k] wcps[i][j, k]]
        # err = KK * v⃗ + b'
        push!(error_vc, vcps[i][:, :] - (-K[1, 1, :, :] .* cys[i][:, :] - K[1, 2, :, :] .* czs[i][:, :]))
        push!(error_wc, wcps[i][:, :] - (-K[2, 1, :, :] .* cys[i][:, :] - K[2, 2, :, :] .* czs[i][:, :]))
    end
    #=
    for i in eachindex(error_vc)
        println("---------")
        print("The average relative error using a local approximation is ")
        println(mean(abs.(error_vc[i] ./ vcps[i])))
        println(" and ")
        println(mean(abs.(error_wc[i] ./ wcps[i])))
        println("-----------")
    end

    for i in eachindex(error_vc)
        println("---------")
        print("The average relative error using a local approximation is ")
        println(mean(abs.(error_vc[i])) ./ mean(abs.(vcps[i])))
        println(" and ")
        println(mean(abs.(error_wc[i])) ./ mean(abs.(wcps[i])))
        println("-----------")
    end
    =#
end

##
# relative differences
mask = zeros(1, 1, size(diffusivities[1])[3:end]...)
# remove dependencies on boundary
for j in 80:(254-80), k in 9:30-9
    mask[1, 1, j, k] = 1.0
end

symmetric_diffusivity = similar(diffusivities)
antisymmetric_diffusivity = similar(diffusivities)
sym_dif = similar(diffusivities[1]) # for temporary allocation
for (i, diffusivity) in enumerate(diffusivities)
    _, _, m, n = size(diffusivity)
    # symmetric part
    for j in 1:m, k in 1:n
        sym_dif[:, :, j, k] .= (diffusivity[:, :, j, k] + diffusivity[:, :, j, k]') * 0.5
    end
    symmetric_diffusivity[i] = copy(sym_dif)

    # anti-symmetric part
    for j in 1:m, k in 1:n
        sym_dif[:, :, j, k] .= (diffusivity[:, :, j, k] - diffusivity[:, :, j, k]') * 0.5
    end
    antisymmetric_diffusivity[i] = copy(sym_dif)
end

error_matrix = zeros(length(diffusivities), length(diffusivities))
sym_distance_matrix = similar(error_matrix)
antisym_distance_matrix = similar(error_matrix)
for i in eachindex(diffusivities), j in eachindex(diffusivities)
    error_matrix[i, j] = norm(diffusivities[i] .* mask - diffusivities[j] .* mask) / norm(0.5 * (diffusivities[i] .* mask + diffusivities[j] .* mask))
    sym_distance_matrix[i, j] = norm(symmetric_diffusivity[i] .* mask - symmetric_diffusivity[j] .* mask) / norm(0.5 * (symmetric_diffusivity[i] .* mask + symmetric_diffusivity[j] .* mask))
    antisym_distance_matrix[i, j] = norm(antisymmetric_diffusivity[i] .* mask - antisymmetric_diffusivity[j] .* mask) / norm(0.5 * (antisymmetric_diffusivity[i] .* mask + antisymmetric_diffusivity[j] .* mask))
end
error_matrix
##

## Apply Spectral Bisection
graph_lap = construct_graph_lap(error_matrix)
weighted_cluster_indices = [collect(1:size(graph_lap)[2])]
weighted_cluster_indices = refine(weighted_cluster_indices, graph_lap)

threshold = heuristic_threshold(error_matrix)
graph_lap = construct_graph_lap(error_matrix, threshold = threshold)
cluster_indices = [collect(1:size(graph_lap)[2])]
cluster_indices = refine(cluster_indices, graph_lap)

##
graph_lap = construct_graph_lap(sym_distance_matrix)
sym_weighted_cluster_indices = [collect(1:size(graph_lap)[2])]
sym_weighted_cluster_indices = refine(sym_weighted_cluster_indices, graph_lap)

threshold = heuristic_threshold(sym_distance_matrix)
graph_lap = construct_graph_lap(sym_distance_matrix, threshold = threshold)
sym_weighted_cluster_indices = [collect(1:size(graph_lap)[2])]
sym_weighted_cluster_indices = refine(sym_weighted_cluster_indices, graph_lap)

##
graph_lap = construct_graph_lap(antisym_distance_matrix)
antisym_weighted_cluster_indices = [collect(1:size(graph_lap)[2])]
antisym_weighted_cluster_indices = refine(antisym_weighted_cluster_indices, graph_lap)


threshold = heuristic_threshold(antisym_distance_matrix)
graph_lap = construct_graph_lap(antisym_distance_matrix, threshold = threshold)
antisym_weighted_cluster_indices = [collect(1:size(graph_lap)[2])]
antisym_weighted_cluster_indices = refine(sym_weighted_cluster_indices, graph_lap)