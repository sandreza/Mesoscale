using JLD2, Plots, Printf, LinearAlgebra
include(pwd() * "/analysis_scripts/" * "post_analysis.jl")

searchdir(path, key) = filter(x -> occursin(key, x), readdir(path))
mesoscale_dir = pwd()
checkpoints = searchdir(mesoscale_dir, "iteration")
tke_list = []
b_list = []
wb_list = []
ub_list = []
vb_list = []
uu_list = []
vv_list = []
ww_list = []

for i in 2:length(checkpoints)
    file = jldopen( mesoscale_dir * "/" * checkpoints[i])
    b = file["tracers"]["b"]["data"][2:end-1, 2:end-1, 2:end-1]
    u = file["velocities"]["u"]["data"][2:end-1, 2:end-1, 2:end-1]
    v = file["velocities"]["v"]["data"][2:end-1, 2:end-1, 2:end-1]
    w = file["velocities"]["w"]["data"][2:end-1, 2:end-1, 2:end-1]
    v = (v[:, 1:end-1, :] + v[:, 2:end, :]) .* 0.5 # put everything on cell centers
    w = (w[:, :, 1:end-1] + w[:,:, 2:end]) .* 0.5  # put everything on cell centers
    close(file)
    n = length(b)
    tke = @. u^2 + v^2 + w^2
    ub = u .* b
    vb = v .* b
    wb = w .* b
    uu = u .* u
    vv = v .* v
    ww = w .* w
    avg_b = sum(b) / n
    avg_tke = sum(tke) / n
    avg_wb = sum(wb) / n
    avg_ub = sum(ub) / n
    avg_vb = sum(vb) / n
    avg_uu = sum(uu) / n
    avg_vv = sum(vv) / n
    avg_ww = sum(ww) / n
    push!(b_list, avg_b)
    push!(tke_list, avg_tke)
    push!(ub_list, avg_ub)
    push!(vb_list, avg_vb)
    push!(wb_list, avg_wb)
    push!(uu_list, avg_uu)
    push!(vv_list, avg_vv)
    push!(ww_list, avg_ww)
    println("the average buoyancy is " * string(avg_b))
    println("the average tke is " * string(avg_tke))
    println("the average wb is " * string(avg_wb))
end

plot(vb_list, label = "vb")
plot!(ub_list, label = "ub")
plot!(wb_list, label = "wb")

##
plot(uu_list[1:end-2], label = "uu")
plot!(vv_list[1:end-2], label = "vv")
plot!(ww_list[1:end-2], label = "ww")