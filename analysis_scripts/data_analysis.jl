using JLD2, Plots, Printf, LinearAlgebra, Statistics
include(pwd() * "/analysis_scripts/" * "post_analysis.jl")
searchdir(path, key) = filter(x -> occursin(key, x), readdir(path))
mesoscale_dir = pwd()
checkpoints = searchdir(mesoscale_dir, "iteration")
println("The checkpoints are", checkpoints)
filename = mesoscale_dir * "/" * checkpoints[end]
println("we are looking at ")
println(filename)
file = jldopen( filename )
b = file["tracers"]["b"]["data"][2:end-1, 2:end-1, 2:end-1]
u = file["velocities"]["u"]["data"][2:end-1, 2:end-1, 2:end-1]
v = file["velocities"]["v"]["data"][2:end-1, 2:end-1, 2:end-1]
w = file["velocities"]["w"]["data"][2:end-1, 2:end-1, 2:end-1]
v = (v[:, 1:end-1, :] + v[:, 2:end, :]) .* 0.5 # put everything on cell centers
w = (w[:, :, 1:end-1] + w[:,:, 2:end]) .* 0.5  # put everything on cell centers
close(file)
n = length(b)
avg_b = sum(b) / n
tke = @. u^2 + v^2 + w^2
avg_tke = sum(tke) / n
println("the average buoyancy is " * string(avg_b))
println("the average tke is " * string(avg_tke))
##
const Nx = size(b)[1]
const Ny = size(b)[2]
const Nz = size(b)[3]
const Lz = 3000.0
const Lx = Ly = 10^6 * 1.0
##
@inline function zC(k)
    return - Lz + (k-0.5) / Nz * Lz
end
@inline function yC(j)
    return  (j-0.5) / Ny * Ly
end
@inline function xC(j)
    return  (j-0.5) / Nx * Lx
end

x = xC.( collect(1:Nx) )
y = yC.( collect(1:Ny) )
z = zC.( collect(1:Nz) )

x2 = (x[1:end-1] + x[2:end]) / 2
y2 = (y[1:end-1] + y[2:end]) / 2
z2 = (z[1:end-1] + z[2:end]) / 2
##
∂x = PartialDerivative((x, 1))
∂y = PartialDerivative((y, 2))
∂z = PartialDerivative((z, 3))
∇ = [∂x, ∂y, ∂z]
##
field_label = [(u, "u"), (v, "v"), (w, "w"), (b, "b"), ( w .* w, "ww")]
selection = 4
field = field_label[selection][1]
label = field_label[selection][2]
cmax = maximum(field)
cmin = minimum(field)
clims = (cmin, cmax)
xind = 96
xloc = @sprintf("%.2f ", x[xind])
p1 = contourf(y, z, field[ xind, :, :]', 
    color = :thermometer, title = "Zonal Slice " * label * " at x=" * xloc,
    xlabel = "Meridional [m]", ylabel = "Depth [m]"
    , clims = clims, linewidth = 0)
##
field_label = [(u, "u"), (v, "v"), (w, "w"), (b, "b"), ( w .* w, "ww"), (v .* b, "vb")]
selection = 5
field = field_label[selection][1]
label = field_label[selection][2]
cmax = maximum(field)
cmin = minimum(field)
clims = (cmin, cmax)
ϕ = sum(field, dims = 1)[1,:,:] ./ Nx
p1 = contourf(y, z, ϕ', 
    color = :thermometer, title = "Zonal Average " * label,
    xlabel = "Meridional [m]", ylabel = "Depth [m]"
    , clims = clims, linewidth = 0)

##
field_label = [(u, "u"), (v, "v"), (w, "w"), (b, "b"), (w .* w, "ww"), (u .* u, "uu"), (v .* v, "vv"), (∂z( w .* w), "d(w .* w) / dz"), (v .* b, "vb"), (∂z(b), "∂z(b)")]
selection = 8
field = sum(field_label[selection][1], dims = (1,2)) ./ (Nx * Ny)
label = field_label[selection][2]
cmax = maximum(field)
cmin = minimum(field)
clims = (cmin, cmax)
p1 = scatter(field[1,1,:],  z,
    color = :blue, xlabel = "Field", label = "Horizontally Averaged " * label,
    ylabel = "Depth [m]", xlims = (cmin, cmax)
    , legend = :bottomright)

##
# day_label = @sprintf("%.2f ", sim_day[i])
# surface values
field_label = [(u, "u"), (v, "v"), (w, "w"), (b, "b"), (u .* b, "ub"), (v .* b , "vb")]
selection = 3 # length(field_label)
field = field_label[selection][1]
label = field_label[selection][2]
cmax = maximum(field)
cmin = minimum(field)
clims = (cmin, cmax)
p1 = contourf(x, y, field[ :, :, end]', 
    color = :thermometer, title = "Surface " * label,
    xlabel = "Zonal [m]", ylabel = "Meridional [m]"
    , clims = clims, linewidth = 0)

## Ertel potential vorticity
U  = [u, v, w]
ω  = ∇ × U
∇b = [∂x(b), ∂y(b), ∂z(b)]
f  = -1e-4
β  = 1e-11
zero_ϕ = ω[1] .* 0
Ω = [zero_ϕ, zero_ϕ, zero_ϕ .+ -1e-4 .+ β .* reshape(y2, (1,191,1))]
total_ω =  ω +  Ω
vec_pv = [total_ω[1] .* ∇b[1] ,  total_ω[2] .* ∇b[2] , total_ω[3] .* ∇b[3]]
pv = vec_pv[1] + vec_pv[2] + vec_pv[3]
magnitude = zeros(size(ω[1]))
@inbounds @simd for i in 1:length(ω[1])
    ω1  = sqrt(total_ω[1][i]^2  + total_ω[2][i]^2  + total_ω[3][i]^2)
    ∇b1 = sqrt(∇b[1][i]^2 + ∇b[2][i]^2 + ∇b[3][i]^2)
    magnitude[i] = ω1 * ∇b1
end
alignment = pv ./ magnitude


##
layer_index = length(z2) - 10
field = ω[3] ./ Ω[3]
label = "instantaneous ω_3 / f "
ϕ = field[ :, :, layer_index]
location_label = @sprintf("%.2f ", z2[layer_index])
cmax = maximum(ϕ)
cmin = minimum(ϕ)
clims = (cmin, cmax)
p1 = contourf(x2, y2, ϕ', 
    color = :thermometer, title = label * " at z=" * location_label * "[m]",
    xlabel = "Zonal [m]", ylabel = "Meridional [m]"
    , clims = clims, linewidth = 0)

field = pv
label = "instantaneous ertel pv "
ϕ = field[ :, :, layer_index]
location_label = @sprintf("%.2f ", z2[layer_index])
cmax = maximum(ϕ)
cmin = minimum(ϕ)
clims = (cmin, cmax)
p2 = contourf(x2, y2, ϕ', 
    color = :thermometer, title = label * " at z=" * location_label * "[m]",
    xlabel = "Zonal [m]", ylabel = "Meridional [m]"
    , clims = clims, linewidth = 0)

field = alignment
label = "instantaneous alignment"
ϕ = field[ :, :, layer_index]
location_label = @sprintf("%.2f ", z2[layer_index])
cmax = maximum(ϕ)
cmin = minimum(ϕ)
clims = (cmin, cmax)
p3 = contourf(x2, y2, ϕ', 
    color = :thermometer, title = label * " at z=" * location_label * "[m]",
    xlabel = "Zonal [m]", ylabel = "Meridional [m]"
    , clims = clims, linewidth = 0)

plot(p1,p2,p3)
##
sf_escale = 32
sf = @. exp( -(y - Ly/2)^2 / (Ly^2 / sf_escale) ) - exp( -(Ly - Ly/2)^2 / (Ly^2 / sf_escale) )
plot(y, sf, ylims = (0,1))

nw_escale = 2*10^4
nw = @. exp( (x - Ly) / nw_escale )
plot!(y, nw, ylims = (0,1))
##

bd_escale = 1 / 200
bd = @. exp(- bd_escale * ( z + Lz))
plot(z, bd, ylims = (0,1))

##
# Checking relaxation profile
Δb = 10 * 2e-3
h = 1000
Nx = 192
Nz = 32
@inline function zC(k)
    return - Lz + (k-0.5) / Nz * Lz
end

@inline function yC(j)
    return  (j-0.5) / Nx * Ly
end

relaxation_profile(j) = Δb * (yC(j)/ Ly)
relaxation_profile_north_2(k) = Δb * ( exp(zC(k)/h) - exp(-Lz/h) ) / (1 - exp(-Lz/h))

tmp1 = relaxation_profile.( collect(1:192))
tmp2 = relaxation_profile_north_2.( collect(1:32))

zonal_average_b = sum(b, dims = (1, )) ./ Nx
tmp1 = reshape(tmp1, (1, 192, 1))
plot( b[1, : , end])
plot!(tmp1[:])
new_field = zonal_average_b .- tmp1
plot(y, new_field[1, :, end])
##
plot( b[2, end , :], legend = false)
for i in 1:10:192
    plot!(b[i, end, :])
end

plot!(tmp2)

##
# Check Geostrophy
∇pʰ = hydrostatic_pressure(b,x,y,z)
zero_ϕ = zeros(size(∇pʰ[1]))
f  = -1e-4
β  = 1e-11
Ω = [zero_ϕ, zero_ϕ, zero_ϕ .+ f .+ β .* reshape(y2, (1,length(y2),1))]
u_avg = (avg_other(u,1)[1:end-1,:,:] + avg_other(u,1)[2:end,:,:]) ./ 2
v_avg = (avg_other(v,1)[1:end-1,:,:] + avg_other(v,1)[2:end,:,:]) ./ 2
w_avg = (avg_other(w,1)[1:end-1,:,:] + avg_other(w,1)[2:end,:,:]) ./ 2
U = [u_avg, v_avg, w_avg]
coriolis_force = [-Ω[3] .* v_avg,  Ω[3] .* u_avg, zero_ϕ]

norm(coriolis_force[1] .+ ∇pʰ[1]) / norm(coriolis_force[1])
norm(coriolis_force[2] .+ ∇pʰ[2]) / norm(coriolis_force[2])

##
field_label = [(coriolis_force[1], "-Ω v"), (-∇pʰ[1], "-∂x(pʰ)"), (coriolis_force[2], "Ω u"), (-∇pʰ[2], "-∂y(pʰ)")]
selection = 2+2*0
field = field_label[selection][1]
label = field_label[selection][2]
cmax = maximum(field)
cmin = minimum(field)
clims = (cmin, cmax)

p1 = contourf(x2, y2, field[ :, :, end]', 
    color = :thermometer, title = "100 [m] depth " * label,
    xlabel = "Zonal [m]", ylabel = "Meridional [m]"
    , clims = clims, linewidth = 0)

selection = 1+2*0
field = field_label[selection][1]
label = field_label[selection][2]
p2 = contourf(x2, y2, field[ :, :, end]', 
    color = :thermometer, title = "100 [m] depth " * label,
    xlabel = "Zonal [m]", yaxis = false
    , clims = clims, linewidth = 0)

field = field_label[selection][1] .- field_label[selection+1][1]
cmax = maximum(field)
cmin = minimum(field)
clims = (cmin, cmax)
label = "deviation"
p3 = contourf(x2, y2, field[ :, :, end]', 
    color = :thermometer, title = "100 [m] depth " * label,
    xlabel = "Zonal [m]", yaxis = false
    , clims = clims, linewidth = 0)
plot(p1,p2,p3)
##
field_label = [(coriolis_force[1], "-Ω v"), (-∇pʰ[1], "-∂x(pʰ)"), (coriolis_force[2], "Ω u"), (-∇pʰ[2], "-∂y(pʰ)")]
selection = 2+2
field = field_label[selection][1]
label = field_label[selection][2]
cmax = maximum(field)
cmin = minimum(field)
clims = (cmin, cmax)

p1 = contourf(x2, y2, field[ :, :, end]', 
    color = :thermometer, title = "100 [m] depth " * label,
    xlabel = "Zonal [m]", ylabel = "Meridional [m]"
    , clims = clims, linewidth = 0)

selection = 1+2
field = field_label[selection][1]
label = field_label[selection][2]
p2 = contourf(x2, y2, field[ :, :, end]', 
    color = :thermometer, title = "100 [m] depth " * label,
    xlabel = "Zonal [m]", yaxis = false
    , clims = clims, linewidth = 0)

field = field_label[selection][1] .- field_label[selection+1][1]
cmax = maximum(field)
cmin = minimum(field)
clims = (cmin, cmax)
label = "deviation"
p3 = contourf(x2, y2, field[ :, :, end]', 
    color = :thermometer, title = "100 [m] depth " * label,
    xlabel = "Zonal [m]", yaxis = false
    , clims = clims, linewidth = 0)
plot(p1,p2,p3)
##
# zonal_average checks 
selection = 2+2
field = mean(field_label[selection][1], dims = 1)[1,:,:]
label = field_label[selection][2]
cmax = maximum(field)
cmin = minimum(field)
clims = (cmin, cmax)

p1 = contourf(y2, z2, field', 
    color = :thermometer, title = "Zonal Average" * label,
    xlabel = "Zonal [m]", ylabel = "Depth [m]"
    , clims = clims, linewidth = 0)

selection = 1+2
field = mean(field_label[selection][1], dims = 1)[1,:,:]
label = field_label[selection][2]
p2 = contourf(y2, z2, field', 
    color = :thermometer, title = "Zonal Average " * label,
    xlabel = "Meriodional [m]", yaxis = false
    , clims = clims, linewidth = 0)

field = mean(field_label[selection][1] .- field_label[selection+1][1], dims = 1)[1,:,:]
# cmax = maximum(field)
# cmin = minimum(field)
clims = (cmin, cmax)
label = "deviation"
p3 = contourf(y2, z2, field', 
    color = :thermometer, title = "Zonal Average " * label,
    xlabel = "Meriodional [m]", yaxis = false
    , clims = clims, linewidth = 0)
plot(p1,p2,p3)
##
selection = 2
field = mean(field_label[selection][1], dims = 1)[1,:,:]
label = field_label[selection][2]
cmax = maximum(field)
cmin = minimum(field)
clims = (cmin, cmax)

p1 = contourf(y2, z2, field', 
    color = :thermometer, title = "Zonal Average" * label,
    xlabel = "Zonal [m]", ylabel = "Depth [m]"
    , clims = clims, linewidth = 0)

selection = 1
field = mean(field_label[selection][1], dims = 1)[1,:,:]
label = field_label[selection][2]
p2 = contourf(y2, z2, field', 
    color = :thermometer, title = "Zonal Average " * label,
    xlabel = "Meriodional [m]", yaxis = false
    , clims = clims, linewidth = 0)

field = mean(field_label[selection][1] .- field_label[selection+1][1], dims = 1)[1,:,:]
# cmax = maximum(field)
# cmin = minimum(field)
clims = (cmin, cmax)
label = "deviation"
p3 = contourf(y2, z2, field', 
    color = :thermometer, title = "Zonal Average " * label,
    xlabel = "Meriodional [m]", yaxis = false
    , clims = clims, linewidth = 0)
plot(p1,p2,p3)
##
gr(size=(700,400))
p3 = plot(p1,p2)
savefig(p3, pwd() * "/geostrophic_0.pdf")
##
p = ∫dz(b, z)
h_p = sum(p, dims = (1,2)) ./ ( (Nx -1) * (Ny-1))
plot(h_p[:], z2)

##
u_prime = u .- mean(u, dims = (1,))
v_prime = v .- mean(v, dims = (1,))
b_prime = b .- mean(b, dims = (1,))
mean_b = mean(b, dims = (1,))
field = b_prime
p1 = contourf(x, y, field[ :, :, end]', 
    color = :thermometer, title = "100 [m] depth ",
    xlabel = "Zonal [m]", ylabel = "Meridional [m]"
    ,  linewidth = 0)

mean_vb = mean(v_prime .* b_prime, dims = (1,))
mean_uv = mean(u_prime .* v_prime, dims = (1,))

p1 = contourf(y, z, mean_vb[1,:,:]', 
    color = :thermometer, title = "<v'b'> ",
    ylabel = "Depth [m]", xlabel = "Meridional [m]"
    ,  linewidth = 0)

    contourf(y, z, mean_vb[1,:,:]', 
    color = :thermometer, title = "<u'v'> ",
    ylabel = "Depth [m]", xlabel = "Meridional [m]"
    ,  linewidth = 0)
##
p2 = contourf(y, z, mean_b[1,:,:]', 
    color = :thermometer, title = "<b> ",
    ylabel = "Depth [m]", xlabel = "Meridional [m]"
    ,  linewidth = 0)
##
∂zb_mean = mean(∂z(b), dims = (1,))
p2 = contourf(y2, z2, ∂zb_mean[1,:,:]', 
    color = :thermometer, title = "<∂z(b)> ",
    ylabel = "Depth [m]", xlabel = "Meridional [m]"
    ,  linewidth = 0)    

𝐟 = reshape(f .+ β * y, (1,192,1))
fvb = 𝐟 .* mean_vb
compare_fvb = avg_other(fvb,1)
##
p3 = contourf(y2, z2, compare_fvb[1,:,:]', 
    color = :thermometer, title = "f<v'b'> ",
    ylabel = "Depth [m]", xlabel = "Meridional [m]"
    ,  linewidth = 0)


plot(p2,p3)

##
ratio = compare_fvb ./ mean(∂zb_mean, dims = (2,)) .* 1027
clims = (minimum(ratio[1,:,21:end]), maximum(ratio[1,:,21:end]))
p3 = contourf(y2, z2, ratio[1,:,:]', 
    color = :thermometer, title = "rho x f x <v'b'>/db/dz ",
    ylabel = "Depth [m]", xlabel = "Meridional [m]"
    ,  linewidth = 0, ylims = (-1000,0), clims = clims)

##
∇b = gradient(b)
b̂ = unitvec(∇b)
norm(sqrt.(b̂[1].^2 .+ b̂[2].^2 .+ b̂[3].^2) .-1)
ba = average_once(b)
ua = average_once(u)
va = average_once(v)
wa = average_once(w)
##
vecu = [ua, va, wa]
vecua = vecu[1] .* b̂[1] + vecu[2] .* b̂[2] + vecu[3] .* b̂[3]
norm(u), norm(v), norm(w)
norm(vecua)
# field = mean(abs.(vecua) , dims=1)[1, :, :]
field = vecua[60,:,:] 
label = " diapycnal velocity "
cmax = maximum(field)
cmin = minimum(field)
clims = (cmin , cmax )

p1 = contourf(y2, z2, field', color = :thermometer, title = "Zonal Slice" * label,
    xlabel = "Zonal [m]", ylabel = "Depth [m]"
    , clims = clims, linewidth = 0)
##
field = vecua[60,:,:] ./ (va[60,:,:] )  .> 0.5
p1 = heatmap(field', color = :blues )
field = vecua[60,:,:] ./ ua[60,:,:] .> 0.5
p2 = heatmap(field', color = :blues )
field = vecua[60,:,:] ./ sqrt.(ua[60,:,:] .^2 .+ va[60,:,:] .^2) .> 0.1
p3 = heatmap(field', color = :blues )
plot(p1,p2,p3)

##
