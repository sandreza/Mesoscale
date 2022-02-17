# find smallest distance to something  per node, Gaurantees that a node has neighbor
# when applying spectral bisection allows us to identify "neighborhoods" in the zero'th eigenmode
function heuristic_threshold(distance_matrix; sort_dim = 2)
    threshold = maximum(sort(distance_matrix, dims = 2)[:, sort_dim]) + eps(maximum(distance_matrix))
    return threshold
end

function construct_graph_lap(distance_matrix; threshold = nothing)
    if isnothing(threshold)
        weighted_graph_lap = Diagonal(sum(distance_matrix, dims = 1)[:]) - distance_matrix
        return weighted_graph_lap
    else
        adj_mat = I - (distance_matrix .< threshold)
        graph_lap = adj_mat - Diagonal(sum(adj_mat, dims = 1)[:])
        return graph_lap
    end
    return nothing
end

# trivial case function, multiple eigenvalue zero, sloppy only need Λ
function trivial_case!(cluster_indices, Λ, V)
    for (i, λ) in enumerate(Λ)
        if abs.(λ) < eps(maximum(Λ))
            indices = eachindex(V[:, i])[abs.(V[:, i]).>eps(maximum(Λ))]
            push!(cluster_indices, indices)
        end
    end
    return nothing
end

# spectral bisection use second eigenvector
function non_trivial_case!(cluster_indices, V)
    bisection_indices = V[:, 2] .> 0
    complement = Bool.(1 .- bisection_indices)
    push!(cluster_indices, eachindex(V[:, 2])[bisection_indices])
    push!(cluster_indices, eachindex(V[:, 2])[complement])
    return nothing
end

# graph indices corresponding to clusters
function graph_clusters!(cluster_indices, graph_lap)
    Λ, V = eigen(graph_lap)
    if sum(Λ .< eps(maximum(Λ))) > 1
        trivial_case!(cluster_indices, Λ, V)
    else
        non_trivial_case!(cluster_indices, V)
    end
end

# given an array of clusters, refine to a further level
function refine(cluster_indices, graph_lap; minimal_group = 2)
    new_cluster_indices = []
    for cluster_index in cluster_indices
        if length(cluster_index) > minimal_group
            new_graph_lap = graph_lap[cluster_index, cluster_index]
            tmp_indices = []
            graph_clusters!(tmp_indices, new_graph_lap)
            for tmp_index in tmp_indices
                push!(new_cluster_indices, cluster_index[tmp_index])
            end
        else
            push!(new_cluster_indices, cluster_index)
        end
    end
    return new_cluster_indices
end
