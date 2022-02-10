hi = 100
if hi % hh == 0
    println("-------")
    println("iteration=", hi)
end
f₀ = 1e-4
inds = collect(1:length(bz[hi, :]))[bz[hi, :].>0]
inds = minimum(inds):maximum(inds)
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