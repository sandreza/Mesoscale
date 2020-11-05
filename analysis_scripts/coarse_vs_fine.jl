using JLD2, Plots, Printf, LinearAlgebra, Statistics
using LaTeXStrings
include(pwd() * "/analysis_scripts/" * "post_analysis.jl")

filename = "/Weno_20_zonal_averages.jld2"
file = jldopen( pwd() * filename)
data_dictionary = get_data(file)
fine_z  = file["grid"]["zC"][2:end-1]
fine_y  = file["grid"]["yC"][2:end-1]
fine_yF = file["grid"]["yF"][2:end-1]
fine_zF = file["grid"]["zF"][2:end-1]
fine_b = mean(data_dictionary["b"][2:end], dims = 1)[1]
fine_u = mean(data_dictionary["u"][2:end], dims = 1)[1]
fine_v = mean(data_dictionary["v"][2:end], dims = 1)[1]
fine_v = (fine_v[2:end,:] + fine_v[1:end-1,:]) ./ 2
data_dictionary = get_data(file)
close(file)

filename = "/Weno_CPU_20_zonal_averages.jld2"
file = jldopen( pwd() * filename)
data_dictionary = get_data(file)
coarse_z  = file["grid"]["zC"][2:end-1]
coarse_y  = file["grid"]["yC"][2:end-1]
coarse_yF = file["grid"]["yF"][2:end-1]
coarse_zF = file["grid"]["zF"][2:end-1]
coarse_b = mean(data_dictionary["b"][2:end], dims = 1)[1]
coarse_u = mean(data_dictionary["u"][2:end], dims = 1)[1]
coarse_v = mean(data_dictionary["v"][2:end], dims = 1)[1]
coarse_v = (coarse_v[2:end,:] + coarse_v[1:end-1,:]) ./ 2
close(file)

filename = "/Weno_CPU_Hack_20_zonal_averages.jld2"
file = jldopen( pwd() * filename)
data_dictionary = get_data(file)
hack_z  = file["grid"]["zC"][2:end-1]
hack_y  = file["grid"]["yC"][2:end-1]
hack_yF = file["grid"]["yF"][2:end-1]
hack_zF = file["grid"]["zF"][2:end-1]
hack_b = mean(data_dictionary["b"][50:end], dims = 1)[1]
hack_u = mean(data_dictionary["u"][50:end], dims = 1)[1]
hack_v = mean(data_dictionary["v"][50:end], dims = 1)[1]
hack_v = (hack_v[2:end,:] + hack_v[1:end-1,:]) ./ 2
close(file)

function plot_field(a, y, yF, z, zF; name = " ")
    py = get_y(a, y, yF) # global scope
    pz = get_z(a, z, zF) # global scope
    clims = (minimum(a), maximum(a))
    further_label = "Zonal and Temporal Average of "
    p1 = contourf(py, pz, a', 
    color = :thermometer, title = " " * name,
    xlabel = "South to North [m]", ylabel = "Depth [m]"
    , clims = clims, linewidth = 0)
    display(p1)
    return p1
end

coarse_plot = plot_field(coarse_b, coarse_y, coarse_yF, coarse_z, coarse_zF, name = "Coarse Buoyancy")
fine_plot = plot_field(fine_b, fine_y, fine_yF, fine_z, fine_zF, name = "Coarse Buoyancy")

function avgy(Φ, n)
    m = size(Φ)[1]
    scale = Int(floor(m/n))
    if ( abs(Int(floor(m/n)) - m/n) > eps(1.0))
        return error
    end
    ny, nz = size(Φ)
    Φ2 = zeros(n, nz)
    for i in 1:n
        Φ2[i,:] .= 0
            for j in 1:scale
                Φ2[i,:] .+= Φ[scale*(i-1) + j,:] / scale
            end
    end
    return Φ2
end

coarse_grained_b = avgy(fine_b, 12)
coarse_grained_plot = plot_field(coarse_grained_b, coarse_y, coarse_yF, coarse_z, coarse_zF, name = "Coarse Grained Buoyancy")

norm(coarse_grained_b - coarse_b) / norm(coarse_grained_b)
norm(coarse_grained_b - hack_b) / norm(coarse_grained_b)

residual_plot = plot_field(coarse_b-coarse_grained_b , coarse_y, coarse_yF, coarse_z, coarse_zF, name = "Slack Buoyancy")