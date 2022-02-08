using JLD2, LinearAlgebra, Statistics

include(pwd() * "/oceananigans_scripts/utils.jl")
include(pwd() * "/diffusivity_scripts/utils_file.jl")


# trial0 and attempt8
diffusivities = []
cases = ["trial1", "trial2", "trial3", "trial4"]

cases = []

for i in 1:5
    push!(cases, "trial"*string(i))
end

for i in 1:7
    push!(cases, "attempt" * string(i))
end

for i in 9:16
    push!(cases, "attempt" * string(i))
end

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

    cs, cys, czs, vcps, wcps, ∇c∇bs, ∇c∇ᵖbs, u⃗c∇bs, u⃗c∇ᵖbs = c_flux_gradient(jl_file, tracer_strings = tracer_strings)
    # to incorporate multiple files, use: [cs[1:2]..., cs[3:4]...]


    tracer_index_partitions = [1:4, 5:8, 9:12]

    for tracer_indices in tracer_index_partitions
        m, n = size(cys[1])

        y, z = get_grid(cs[1], jl_file; ghost = 3)
        y̅ = avg(y)
        z̅ = avg(z)

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
            for i in tracer_indices
                ωⁱ = case_weights[i]
                # calculate matrix entries
                A[1, 1] += ωⁱ * cys[i][j, k] * cys[i][j, k]
                A[1, 2] += ωⁱ * cys[i][j, k] * czs[i][j, k]
                A[2, 2] += ωⁱ * czs[i][j, k] * czs[i][j, k]
                # calculate rhs 1 for κ¹¹, κ¹²
                b1[1] -= ωⁱ * cys[i][j, k] * vcps[i][j, k]
                b1[2] -= ωⁱ * czs[i][j, k] * vcps[i][j, k]
                # calculate rhs 2 for κ²¹, κ²²
                b2[1] -= ωⁱ * cys[i][j, k] * wcps[i][j, k]
                b2[2] -= ωⁱ * czs[i][j, k] * wcps[i][j, k]
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
            # println("finishing j = $j and k = $k")
        end

        # Compute Error on a Per Tracer Basis
        error_vc = []
        error_wc = []
        for i in tracer_indices
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
        push!(diffusivities, K)
    end
end
##
# relative differences
error_matrix = zeros(length(diffusivities), length(diffusivities))
for i in eachindex(diffusivities), j in eachindex(diffusivities)
    error_matrix[i,j] = norm(diffusivities[i] - diffusivities[j]) / mean(norm.(diffusivities))
    # println(norm(diffusivities[i] - diffusivities[j]) / mean(norm.(diffusivities)))
end
error_matrix