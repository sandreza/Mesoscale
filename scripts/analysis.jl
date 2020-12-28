using JLD2, LinearAlgebra, Oceananigans
include(pwd() * "/scripts/vizinanigans.jl")
include(pwd() * "/analysis_scripts/" * "post_analysis.jl") # Gradients etc. here
filename = pwd() * "/Channel_16_checkpoint_iteration6317902.jld2"
file = jldopen( filename )
grid = file["grid"]
xS, yS, zS = size(grid)
btmp = file["tracers"]["b"]["data"]
ghost = Int((length(btmp[:,1,1])-xS)/2)
b = file["tracers"]["b"]["data"][ghost:(xS + ghost - 1), ghost:(yS + ghost - 1), ghost:(zS + ghost - 1)]
u = file["velocities"]["u"]["data"][ghost:(xS + ghost - 1), ghost:(yS + ghost - 1), ghost:(zS + ghost - 1)]
v = file["velocities"]["v"]["data"][ghost:(xS + ghost - 1), ghost:(yS + 1 + ghost - 1), ghost:(zS + ghost - 1)]
w = file["velocities"]["w"]["data"][ghost:(xS + ghost - 1), ghost:(yS + ghost - 1), ghost:(zS + 1 + ghost - 1)]

x = grid.xC[1:xS]
y = grid.yC[1:yS]
z = grid.zC[1:zS]
# Cell-Centered Values
v = (v[:, 1:end-1, :] + v[:, 2:end, :]) .* 0.5
w = (w[:, :, 1:end-1] + w[:,:, 2:end]) .* 0.5
close(file)
## Create Gradient
∂x = PartialDerivative((x, 1))
∂y = PartialDerivative((y, 2))
∂z = PartialDerivative((z, 3))
∇ = [∂x, ∂y, ∂z]
ω = ∇ × [u, v, w]
∇b = [∂x(b), ∂y(b), ∂z(b)]
coriolis =  -1e-4 .+ 1e-11 * reshape(y, (1, 192 ,1)) 
fpβ = coriolis[1:end-1] + coriolis[2:end] 
ω∇b = ω⋅∇b
pv = ω∇b + fpβ .* ∇b[3]
alignment = ω∇b ./ ( sqrt.(ω⋅ω) .* sqrt.(∇b⋅∇b) )
states = [u, abs.(u), v, w, b, ω[3], ω[1], ω[2], abs.(ω∇b), abs.(pv)]
statenames=  ["u", "|u|", "v", "w", "b", "ω₃", "ω₁", "ω₂", "|ω⋅∇b|",  "|Ertel PV|"]
visualize(states, statenames = statenames, quantiles = (0.10, 0.99))

##