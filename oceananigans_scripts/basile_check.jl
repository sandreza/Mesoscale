using JLD2
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

# eigen(Symmetric())
f₀ = 1e-4
N² = bz[hi, zi_b:zi_s]
mat = sturm_louiville_operator(f₀^2 ./ N², Δz)
λ⁻², V = eigen(mat)
λ₁ = (-λ⁻²[end-1])^(-0.5) # First Rossby Deformation Radius
V₁ = V[:, end-1] # first baroclinic mode
λ₂ = (-λ⁻²[end-2])^(-0.5) # Second Rossby Deformation Radius
V₁ = V[:, end-1] # second baroclinic mode

P = V' # projection operator since orthonormal
B₁ = P[end-1, :]' # last row is projection operator

uu = u[hi, zi_b-1:zi_s]
û = P * uu
norm(û[end] * V[:, end] - uu) / norm(uu)
norm(û[end] * V[:, end] + û[end-1] * V[:, end-1] - uu) / norm(uu)
norm(û[end] * V[:, end] + û[end-1] * V[:, end-1] + û[end-2] * V[:, end-2] - uu) / norm(uu)
û = reverse(abs.(P * uu[hi, zi_b-1:zi_s]))

