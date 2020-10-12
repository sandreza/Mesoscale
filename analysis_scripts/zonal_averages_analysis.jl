using JLD2, Plots, Printf, LinearAlgebra, Statistics
using LaTeXStrings
include(pwd() * "/analysis_scripts/" * "post_analysis.jl")


filename = "/Abernathy_20_zonal_averages.jld2"
file = jldopen( pwd() * filename)
z = file["grid"]["zC"][2:end-1]
y = file["grid"]["yC"][2:end-1]
yF = file["grid"]["yF"][2:end-1]
zF = file["grid"]["zF"][2:end-1]
data_dictionary = get_data(file)
close(file)
##
function get_y(a, y, yF)
    length(a[:,1]) == length(y) ? y : yF
end
function get_z(a, z, zF)
    length(a[1,:]) == length(z) ? z : zF
end

function plot_field(a; name = " ")
    py = get_y(a, y, yF) # global scope
    pz = get_z(a, z, zF) # global scope
    clims = (minimum(a), maximum(a))
    p1 = contourf(py, pz, a', 
    color = :thermometer, title = "Zonal and Temporal Average of " * name,
    xlabel = "Meridional [m]", ylabel = "Depth [m]"
    , clims = clims, linewidth = 0)
    display(p1)
end
##
tmp = norm(data_dictionary["u"][end] - data_dictionary["u"][end-1])
tmp /= norm(data_dictionary["u"][end])

vb = mean(data_dictionary["vb"][end-20:end])
vb = (vb[2:end, :] + vb[1:end-1, :]) ./ 2
b = mean(data_dictionary["b"][end-20:end])
# plot_field(b, name = "b")
u = mean(data_dictionary["u"][end-20:end])
v = mean(data_dictionary["v"][end-20:end])
v = (v[2:end, :] + v[1:end-1, :]) ./ 2
vpbp = vb - v .* b
f  = -1e-4
β  = 1e-11
𝐟 = reshape(f .+ β * y, (192,1))
ρ = 1027
top = ρ .* 𝐟 .* vpbp
Δz = reshape(z[2:end] - z[1:end-1], (1, length(z)-1))
b = (b[:,2:end] + b[:,1:end-1]) ./ 2

bottom = mean(b ./ Δz, dims = 1)
top = (top[:,1:end-1] + top[:, 2:end]) ./ 2

##
py = y
pz = (z[1:end-1] + z[2:end]) ./ 2
a = top ./ bottom
a = a[:, 16:end]
pz = pz[16:end]
clims = (minimum(a), maximum(a))
p1 = contourf(py, pz, a', 
    color = :thermometer, title = "Zonal and Temporal Average " * " rho x f x <v'b'>/db/dz",
    xlabel = "Meridional [m]", ylabel = "Depth [m]"
    , clims = clims, linewidth = 0)



##
b = mean(data_dictionary["b"][2:end])
plot_field(b, name = "b")

##

py = y
pz = (z[1:end-1] + z[2:end]) ./ 2
a = top ./ bottom
a = mean(a[:, 16:end], dims = 1)
pz = pz[16:end]
clims = (minimum(a), maximum(a))
p1 = plot(a[:], pz)
plot(bottom[16:end], pz)

##
u = mean(data_dictionary["u"][2:end])

plot_field(u)

u_avg = mean(u[20:end-20, :], dims = 1)

p1 = plot(u_avg[:], z, title = "horizontal average", 
label = L"$ \overline{u} $", xlabel = "speed [m/s]", ylabel = "depth",
legend = :bottomright)
##

label = "b"
v = mean(data_dictionary[label][2:end])
v_avg = mean(v[20:end-20, :], dims = 1)
p2 = plot(v_avg[:], z, title = "horizontal average", 
label = L"$\overline{b}$", xlabel = "bouyancy [m/s²]", ylabel = "depth",
legend = :bottomright, color = :red)


v_avg = mean(vpbp[20:end-20, :], dims = 1)
p4 = plot(v_avg[:] *10^4, z, title = "horizontal average", 
label = L"$\overline{v' b'}$", xlabel = "[cm²/s³]", ylabel = "depth",
legend = :bottomleft, color = :orange)

##

a = top ./ bottom
z2 = (z[2:end] + z[1:end-1]) / 2
v_avg = mean(a, dims = 1)
p3 = plot(v_avg[:], z2, title = "horizontal average", 
label = L"$ \rho  \overline{f v'b'} / \overline{\partial_z b} $", xlabel = "stress [N/m]", ylabel = "depth",
legend = :right, color = :purple)


p8 = plot(p1,p2,p3,p4)
savefig(p8, "plots.pdf")


