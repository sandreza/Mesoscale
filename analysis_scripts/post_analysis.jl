function ∂(field, grid, direction)
    gp = length(grid)
    access = define_shifts(1:gp-1, :, direction)
    access_shift = define_shifts(2:gp, :, direction)
    Δf = field[access...] - field[access_shift...]
    Δg = grid[1:end-1] - grid[2:end]
    reshape_Δg = define_shifts(gp-1, 1, direction)
    Δg = reshape(Δg, reshape_Δg)
    return Δf ./ Δg 
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
