function ∂(field, grid, direction)
    gp = length(grid)
    access = define_shifts(1:gp-1, :, direction)
    access_shift = define_shifts(2:gp, :, direction)
    Δf = field[access_shift...] - field[access...] 
    Δg = grid[2:end] - grid[1:end-1]
    reshape_Δg = define_shifts(gp-1, 1, direction)
    Δg = reshape(Δg, reshape_Δg)
    return Δf ./ Δg 
end

function ∫dz(field, grid; direction = 3)
    gp = length(grid)
    access = define_shifts(1:gp-1, :, direction)
    access_shift = define_shifts(2:gp, :, direction)
    f_avg = (field[access...] + field[access_shift...]) ./ 2
    Δg = grid[2:end] - grid[1:end-1]
    reshape_Δg = define_shifts(gp-1, 1, direction)
    Δg = reshape(Δg, reshape_Δg)
    ∫f = zeros(size(field[access...]))
    ∫f[:, :, 1] = Δg[1] * field[:, :, 1]
    for i in 2:length(Δg)
        ∫f[:, :, i] = Δg[i] * field[:, :, i] + ∫f[:, :, i-1]
    end
    return ∫f
end

function define_shifts(tmp, tmp2, direction)
    if direction==1
        return (tmp, tmp2, tmp2)
    elseif direction==2
        return (tmp2, tmp, tmp2)
    else
        return (tmp2, tmp2, tmp)
    end
    return nothing
end
#=
# this instead?
function appropriate_dims(n1, n2, N)
    y = ones(Int, n1)
    y[n2] = N
    return Tuple(y)
end
=# 
function avg_other(a, direction)
    nx, ny, nz = size(a)
    if direction==1
        b = (a[:, 1:ny-1, :] + a[:, 2:ny, :]) .* 0.5
        c = (b[:, :, 1:nz-1] + b[:, :, 2:nz]) .* 0.5
        return c
    elseif direction==2    
        b = (a[:, :, 1:nz-1] +a[:, :, 2:nz])  .* 0.5
        c = (b[1:nx-1, :, :] + b[2:nx, :, :]) .* 0.5
        return c
    else
        b = (a[1:nx-1, :, :] + a[2:nx, :, :]) .* 0.5
        c = (b[:, 1:ny-1, :] + b[:, 2:ny, :]) .* 0.5
        return c
    end
    return nothing
end



import Base: *

struct PartialDerivative{T}
    d::T
end
struct Gradient{T}
    g::T
end
struct Field{T}
    data::T
end

function *(∇::PartialDerivative, ϕ::AbstractArray)
    tmp = ∂(ϕ, ∇.d[1], ∇.d[2])
    tmp2 = avg_other(tmp, ∇.d[2])
    return tmp2
end

function (p::PartialDerivative)(ϕ::AbstractArray)
    return *(p, ϕ)
end

*(∇::PartialDerivative, ϕ::Field) = *(∇, ϕ.data)

function hydrostatic_pressure(b,x,y,z)
    p = ∫dz(b, z)
    ∂ᶻp = (b[:,:,1:end-1] + b[:,:,2:end]) ./ 2 
    ∂ᶻp = avg_other(∂ᶻp, 3)
    Δx = reshape(x[1:end-1] - x[2:end], (length(x)-1, 1, 1))
    ∂ˣp = (p[1:end-1,:,:] - p[2:end,:,:]) ./ Δx
    ∂ˣp = (∂ˣp[:,1:end-1,:] + ∂ˣp[:,2:end,:]) ./ 2
    Δy = reshape(x[1:end-1] - x[2:end], (1,length(y)-1, 1))
    ∂ʸp = (p[:,1:end-1,:] - p[:,2:end,:]) ./ Δy
    ∂ʸp = (∂ʸp[1:end-1,:,:] + ∂ʸp[2:end,:,:]) ./ 2
    ∇p = [∂ˣp, ∂ʸp, ∂ᶻp]
    return ∇p
end

import LinearAlgebra: dot
function dot(ω::Array{Array{Float64,3},1}, ∇b::Array{Array{Float64,3},1}) 
    tmparray = copy(ω[1]) .* 0.0
    for i in eachindex(ω)
        tmparray += ω[i] .* ∇b[i]
    end
    return tmparray
end
