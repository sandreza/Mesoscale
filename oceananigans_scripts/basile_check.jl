using JLD2, Statistics
using LinearAlgebra
include("utils.jl")
jl_file = jldopen("relaxation_channel_tracers_restarted_smooth_forcing_case_attempt0_averages.jld2", "r+")

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

vbp = avgy(avgz(avgy(vb) - avgy(v) .* b))
wbp = avgy(avgz(avgz(wb) - avgz(w) .* b))

u = get_field("u", time_index, jl_file)
v = get_field("v", time_index, jl_file)
w = get_field("w", time_index, jl_file)

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

𝒟(c₁, c₂) = f₀ * c₁ * 𝒰 * λ₁ * exp(c₂ * 𝒰 / λ₁)

𝒰list = []
λlist = []
D̅list = []
Dlist = []
δD̅list = []
wkb_λlist = []
V₁list = []

for hi in 1:255
    hh = 10
    if hi % hh == 0
        println("-------")
        println("iteration=", hi)
    end
    f₀ = 1e-4
    inds = collect(1:length(bz[hi, :]))[bz[hi, :].>0]
    N² = bz[hi, :] # just take all slopes

    mat = sturm_louiville_operator(f₀^2 ./ N², Δz)
    λ⁻², V = eigen(mat)
    λ⁻² = real.(λ⁻²) # should be positive definite anyway
    λ₀ind = argmin(abs.(λ⁻²))
    λ₁ = (-λ⁻²[λ₀ind-1])^(-0.5) # First Rossby DeformatiV₁ = V[:, λ₀ind-1] # first baroclinic mode
    λ₂ = (-λ⁻²[λ₀ind-2])^(-0.5) # Second Rossby Deformation Radius

    V₀ = V[:, λ₀ind] # barotropic mode
    V₁ = V[:, λ₀ind-1] # first baroclinic mode
    V₂ = V[:, λ₀ind-2] # second baroclinic mode
    wkb_λ₁ = mean(sqrt.(abs.(N²))) / (π * f₀) * 3000

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
    if hi % hh == 0
        println("er1=", er1)
        println("er2=", er2)
        println("er3=", er3)
        println("deformation ", λ₁)
        println("deformation wkb ", wkb_λ₁)
        println("-------")
    end

    𝒰 = abs(û[end-1])

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
end

##
𝒟(c₁, c₂; 𝒰 = 𝒰, λ₁ = λ₁, f₀ = -1e-4) = f₀ * c₁ * 𝒰 * λ₁ * exp(c₂ * 𝒰 / λ₁)
vpbp = avgz(avgy(vb) - avgy(v) .* b)

# c₁ is scaled by 100
# c₂ is scaled by hours

function loss_function(c₁, c₂; norm_function = maximum, zinds = 10:28, yinds = 80:180)
    loss = []

    for hi in yinds
        𝒰 = 𝒰list[hi]
        λ₁ = λlist[hi]
        push!(loss, norm(𝒟(c₁, 86400 / 24 * c₂, 𝒰 = 𝒰, λ₁ = λ₁) * uz[hi, zinds] - vpbp[hi, zinds]))
    end
    return norm_function(loss)
end

NN = 100
c1vals = reshape(collect(0:NN) ./ NN * 1.0, (NN + 1, 1))
c2vals = reshape(collect(0:NN) ./ NN * 1.0, (1, NN + 1))
lossvals = loss_function.(c1vals, c2vals)

println("relative decrease ", minimum(lossvals) / loss_function(0, 0))

zinds = 10:28
hi = 100
norm(𝒟(0.5, 0.03, 𝒰 = 𝒰, λ₁ = λ₁) * uz[hi, zinds] - vpbp[hi, zinds]) / norm(vpbp[hi, zinds])

c1min = c1vals[argmin(lossvals)[1]]
c2min = c2vals[argmin(lossvals)[2]]

for hi in 70:5:180
    println("---------")
    println(" y = ", y[hi])
    er1 = norm(𝒟(c1min, c2min, 𝒰 = 𝒰, λ₁ = λ₁) * uz[hi, zinds] - vpbp[hi, zinds]) / norm(vpbp[hi, zinds])
    println("error = ", er1)
    println("---------")
end
