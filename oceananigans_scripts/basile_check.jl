using JLD2, Statistics
using LinearAlgebra
include("utils.jl")
jl_file = jldopen("relaxation_channel_tracers_restarted_smooth_forcing_case_attempt0_averages.jld2", "r+")
jl_file = jldopen("fluxernathy_tracers_restarted_smooth_forcing_averages.jld2", "r+")
##
keys(jl_file)
keys(jl_file["timeseries"])
time_keys = keys(jl_file["timeseries"]["vb"])

##
time_index = 22
b = get_field("b", time_index, jl_file)
y, z = get_grid(b, jl_file; ghost = 3)

y̅ = avg(y)
z̅ = avg(z)
bz_p = ∂z(b, z)
by_p = ∂y(b, y)

bz = avgy(bz_p)
by = avgz(by_p)

vb = get_field("vb", time_index, jl_file)
wb = get_field("wb", time_index, jl_file)
w = get_field("w", time_index, jl_file)

u = get_field("u", time_index, jl_file)
v = get_field("v", time_index, jl_file)
w = get_field("w", time_index, jl_file)

vbp = avgy(avgz(avgy(vb) - avgy(v) .* b))
wbp = avgy(avgz(avgz(wb) - avgz(w) .* b))

uz = ∂z(u, z)

##
hi = 100 # horizontal index
# large zi means surface
zi_s = 31  # surface index end
zi_b = 5   # bottom index start
println("looking at location y= ", y̅[hi], " from depth ", " z = ", z̅[zi_b], " to z = ", z̅[zi_s])
bz[hi, zi_b:zi_s]

## 
# Step 1: Figure out how much energy is in each of the modes
# only use indices 2 onwards of κ, thus defined at cell faces
function sturm_louiville_operator(κ, Δz)
    N = length(κ)
    matrix = Tridiagonal(zeros(N + 1, N + 1))
    for i in 2:N
        matrix.d[i] = -2 * (κ[i-1] + κ[i]) * 0.5
        matrix.dl[i] = 1 * κ[i]
        matrix.du[i] = 1 * κ[i]
    end
    matrix.d[1] = -κ[1]
    matrix.dl[1] = κ[1]
    matrix.d[end] = -κ[N]
    matrix.dl[end] = κ[N]
    matrix *= (Δz)^(-2)
    matrix.du .= matrix.dl
    return Symmetric(Array(matrix))
end

##
Δz = z[2] - z[1]
mat = sturm_louiville_operator(collect(1:10), Δz)
eigen(mat)

𝒰list = Float64[]
practicalU_list = Float64[]
λlist = Float64[]
D̅list = []
Dlist = []
δD̅list = []
wkb_λlist = Float64[]
V₁list = []

for hi in 1:255
    hh = 10
    if hi % hh == 0
        println("-------")
        println("iteration=", hi)
    end
    f₀ = -1e-4 + 1e-11 * y[hi]
    H = 3000
    inds = collect(1:length(bz[hi, :]))[bz[hi, :].>0]
    N² = bz[hi, :] # just take all slopes

    mat = sturm_louiville_operator(f₀^2 ./ N², Δz)
    λ⁻², V = eigen(mat)
    λ⁻² = real.(λ⁻²) # should be positive definite anyway
    λ₀ind = argmin(abs.(λ⁻²))
    λ₁ = (-λ⁻²[λ₀ind-1])^(-0.5) # First Rossby Deformation Radius V₁ = V[:, λ₀ind-1] # first baroclinic mode
    λ₂ = (-λ⁻²[λ₀ind-2])^(-0.5) # Second Rossby Deformation Radius

    V₀ = V[:, λ₀ind] # barotropic mode
    V₁ = V[:, λ₀ind-1] # first baroclinic mode
    V₂ = V[:, λ₀ind-2] # second baroclinic mode
    wkb_λ₁ = mean(sqrt.(abs.(N²))) / (π * abs(f₀)) * H

    P = inv(V) # projection operator 
    B₁ = P[λ₀ind-1, :]' # last row is projection operator

    # extract velocity component
    uu = u[hi, :]

    # project
    û = P * uu
    er1 = norm(û[λ₀ind] * V₀ - uu) / norm(uu)
    er2 = norm(û[λ₀ind] * V₀ + û[λ₀ind-1] * V₁ - uu) / norm(uu)
    er3 = norm(û[λ₀ind] * V₀ + û[λ₀ind-1] * V₁ + û[λ₀ind-2] * V₂ - uu) / norm(uu)
    # û = reverse(abs.(P * uu))

    𝒰 = abs(û[λ₀ind-1])
    Upractical = norm(u[hi, :] .- mean(u[hi, :]))

    if hi % hh == 0
        println("er1=", er1)
        println("er2=", er2)
        println("er3=", er3)
        println("deformation ", λ₁)
        println("deformation wkb ", wkb_λ₁)
        println("𝒰 = ", 𝒰)
        println("Quick Estimate = ", Upractical)
        println("-------")
    end



    f₀ = 1e-4


    avgz(vb)[hi, :] ./ uz[hi, :]

    vbp[hi, :]

    D = vbp[hi, :] ./ uz[hi, :]
    D̅ = mean(-D[D.<0])
    δD̅ = std(-D[D.<0])

    push!(𝒰list, 𝒰)
    push!(λlist, λ₁)
    push!(D̅list, D̅)
    push!(δD̅list, δD̅)
    push!(Dlist, D)
    push!(wkb_λlist, wkb_λ₁)
    push!(V₁list, V₁)
    push!(practicalU_list, Upractical)
end

##
𝒟(c₁, c₂; 𝒰 = 𝒰, λ₁ = λ₁, f₀ = -1e-4) = sign(f₀) * c₁ * 𝒰 * λ₁ * exp(c₂ * 𝒰 / λ₁)
vpbp = avgz(avgy(vb) - avgy(v) .* b)

# c₁ is scaled by 100
# c₂ is scaled by hours

function loss_function(c₁, c₂; norm_function = maximum, zinds = 10:28, yinds = 80:180, debug = false)
    loss = []

    for hi in yinds
        𝒰 = abs(𝒰list[hi])
        λ₁ = λlist[hi]
        f₀ = -1e-4 + 1e-11 * y[hi]
        push!(loss, norm(𝒟(c₁, c₂, 𝒰 = 𝒰, λ₁ = λ₁) * uz[hi, zinds] * f₀ + vpbp[hi, zinds]))
    end
    if debug
        return loss
    else
        return norm_function(loss)
    end
end

NN = 100
# A priori determine ranges for optimization: 
zinds = 10:28
hi = 100 # middle horizontal index
𝒰1 = 𝒰list[hi]
λ₁1 = λlist[hi]
f₀ = -1e-4 + 1e-11 * y[hi]
p¹ = mean(-vpbp[hi, zinds] ./ (f₀ * uz[zinds]))
hi = 150
𝒰2 = 𝒰list[hi]
λ₁2 = λlist[hi]
f₀ = -1e-4 + 1e-11 * y[hi]
p² = mean(-vpbp[hi, zinds] ./ (f₀ * uz[hi, zinds]))
# p¹ =  c₁ * 𝒰1 * λ₁1 * exp(c₂ * 𝒰1 / λ₁1)
# p² =  c₁ * 𝒰2 * λ₁2 * exp(c₂ * 𝒰2 / λ₁2)
# c₂ = ln(p¹ * 𝒰2 * λ₁2 / (p²) / (𝒰1 / λ₁1 -  𝒰2 / λ₁2)
c₂guess = log((p¹ * 𝒰2 * λ₁2) / (p² * 𝒰1 * λ₁1)) / (𝒰1 / λ₁1 - 𝒰2 / λ₁2)
c₁guess = -p¹ / (𝒰1 * λ₁1 * exp(c₂guess * 𝒰1 / λ₁1))
-p² / (𝒰2 * λ₁2 * exp(c₂guess * 𝒰2 / λ₁2))

c1vals = reshape(collect(0:NN) ./ NN * 5 * c₁guess, (NN + 1, 1))
c2vals = reshape(collect(0:NN) ./ NN * 5 * c₂guess, (1, NN + 1))
lossvals = loss_function.(c1vals, c2vals)
loss_function(c₁guess, c₂guess, debug = true)
loss_function(0, 0)
println("relative decrease ", minimum(lossvals) / loss_function(0, 0))

zinds = 10:28
hi = 100
𝒰 = 𝒰list[hi]
λ₁ = λlist[hi]
f₀ = -1e-4 + 1e-11 * y[hi]
norm(𝒟(c₁guess, c₂guess, 𝒰 = 𝒰, λ₁ = λ₁) * uz[hi, zinds] * f₀ + vpbp[hi, zinds]) / norm(vpbp[hi, zinds])

c1min = c1vals[argmin(lossvals)[1]]
c2min = c2vals[argmin(lossvals)[2]]
norm(𝒟(c1min, c2min, 𝒰 = 𝒰, λ₁ = λ₁) * uz[hi, zinds] - vpbp[hi, zinds]) / norm(vpbp[hi, zinds])
##
# 0.2 
# c2min = 4.0e5 # 2.4e5 
for hi in 70:5:180
    println("---------")
    println(" y = ", y[hi])
    println(" distance percent ", y[hi] / y[end])
    𝒰 = practicalU_list[hi] # 𝒰list[hi]
    λ₁ = wkb_λlist[hi] # λlist[hi]
    f₀ = -1e-4 + 1e-11 * y[hi]
    er1 = norm(𝒟(c1min, c2min, 𝒰 = 𝒰, λ₁ = λ₁) * uz[hi, zinds] * f₀ + vpbp[hi, zinds]) / norm(vpbp[hi, zinds])
    println("error = ", er1)
    println("---------")
end

##
dlist = Float64[]
truedlist = Float64[]
trued2list = Float64[]

avguz = avgy(uz)

for hi in 1:255
    println("---------")
    f₀ = -1e-4 + 1e-11 * y[hi]
    𝒰 = practicalU_list[hi] # 𝒰list[hi]
    λ₁ = wkb_λlist[hi] # λlist[hi]
    push!(dlist, -𝒟(c1min, c2min, 𝒰 = 𝒰, λ₁ = λ₁))
    push!(truedlist, mean(vpbp[hi, zinds] ./ uz[hi, zinds]) ./ f₀)
    push!(trued2list, mean(-vpbp[hi, zinds] ./ by[hi, zinds]))
    println("Thermal wind error = ", norm(by[hi, zinds] ./ avguz[hi, zinds] .- (-f₀), Inf) / abs(f₀))
end

options = (; ylabel = "K", xlabel = "y [km]", ylabelsize = 32,
    xlabelsize = 32, xgridstyle = :dash, ygridstyle = :dash, xtickalign = 1,
    xticksize = 30, ytickalign = 1, yticksize = 30,
    xticklabelsize = 30, yticklabelsize = 30)

fig = Figure(resolution = (1800, 1300), title = "Diffusivity Comparison")
titlestring = "Diffusivity Comparison"
ax1 = Axis(fig[1, 1]; options..., title = titlestring, titlesize = 30)

plot_string_1 = "Simulation Estimate 1"
plot_string_2 = "Simulation Estimate 2"
plot_string_3 = "c₁ = 0.2 [*], c₂ = 2.4e5 [s]"

ln1 = lines!(ax1, y̅ ./ 1e3, truedlist, color = :black)
ln2 = lines!(ax1, y̅ ./ 1e3, trued2list, color = :purple)
ln3 = lines!(ax1, y̅ ./ 1e3, dlist, color = :red)
axislegend(ax1, [ln1, ln2, ln3], [plot_string_1, plot_string_2, plot_string_3], position = :rt, labelsize = 30)


##
options = (; xlabel = "y [km]", ylabelsize = 32,
    xlabelsize = 32, xgridstyle = :dash, ygridstyle = :dash, xtickalign = 1,
    xticksize = 30, ytickalign = 1, yticksize = 30,
    xticklabelsize = 30, yticklabelsize = 30)

fig2 = Figure(resolution = (1800, 1300))
titlestring_1 = "Rossby Deformation Radius"
ax1 = Axis(fig2[1, 1]; options..., title = titlestring_1, titlesize = 30, ylabel = "λ [m]")

titlestring_2 = "Baroclinic Mode Amplitue"
ax2 = Axis(fig2[2, 1]; options..., title = titlestring_2, titlesize = 30, ylabel = "𝒰 [m/s]")

plot_string_1 = "Eigenvalue Solve"
plot_string_2 = "WKB Estimate"

ln1 = lines!(ax1, y̅ ./ 1e3, λlist, color = :black)
ln2 = lines!(ax1, y̅ ./ 1e3, wkb_λlist, color = :red)

axislegend(ax1, [ln1, ln2], [plot_string_1, plot_string_2], position = :rb, labelsize = 30)

plot_string_1 = "Sturm-Liouville Solve + Projection"
plot_string_2 = "L² Norm of Velocity Perturbation"

ln1 = lines!(ax2, y̅ ./ 1e3, 𝒰list, color = :black)
ln2 = lines!(ax2, y̅ ./ 1e3, practicalU_list, color = :red)

axislegend(ax2, [ln1, ln2], [plot_string_1, plot_string_2], position = :cb, labelsize = 30)

