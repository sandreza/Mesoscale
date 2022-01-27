using LinearAlgebra 
using SparseArrays
using BenchmarkTools

N = 32
d = 2 * ones(N)
dl = -ones(N-1)
TM = Tridiagonal(dl, d, dl)
i3 = TM .* 0 + I
TM[1] = TM[end] = 1

Δ = kron(i3, TM) + kron(TM, i3) + I 
s_Δ = sparse(Δ)
ϕ = zeros(N,N)

ϕ .= 0.0
ϕ[:,1] .= 1
φ = copy(ϕ )

s_Δ = sparse(Δ)

cΔ =  cholesky(Δ)
c_Δ =  cholesky(s_Δ)
answ1 = Δ \ ϕ[:]
answ2 = c_Δ \ ϕ[:]
norm(answ1 - answ2)
@benchmark c_Δ \ ϕ[:]
