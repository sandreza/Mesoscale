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
    Δy = reshape(y[1:end-1] - y[2:end], (1,length(y)-1, 1))
    ∂ʸp = (p[:,1:end-1,:] - p[:,2:end,:]) ./ Δy
    ∂ʸp = (∂ʸp[1:end-1,:,:] + ∂ʸp[2:end,:,:]) ./ 2
    ∇p = [∂ˣp, ∂ʸp, ∂ᶻp]
    return ∇p
end

function get_data(file)
    data_dictionary = Dict()
    for key in keys(file["timeseries"])
        println(key)
        tmp = []
        for t in keys(file["timeseries"]["t"])
            if key != "t"
                push!(tmp, file["timeseries"][key][t][1, :,:])
            else
                push!(tmp, file["timeseries"][key][t])
            end
        end
        data_dictionary[key] = tmp
    end
    return data_dictionary
end


##
# Get diapycnal component
function gradient(field)
    return [∂x(field), ∂y(field), ∂z(field)]
end

function unitvec(∇b)
    ∂ˣb, ∂ʸb, ∂ᶻb = ∇b
    magnitude = copy(∂ˣb)
    for i in 1:length(∂ˣb)
        magnitude[i] = sqrt(∂ˣb[i]^2 + ∂ʸb[i]^2 + ∂ᶻb[i]^2)
    end
    return [∂ˣb ./ magnitude, ∂ʸb ./ magnitude, ∂ᶻb ./ magnitude]
end

function average_once(b)
    c = (b[2:end, :, :] + b[1:end-1, :, :]) ./2
    c = (c[:, 2:end, :] + c[:, 1:end-1 , :]) ./2
    c = (c[:, :, 2:end] + c[:, :, 1:end-1]) ./2
    return c
end
    
## Coarse-grained
function avg(Φ, n)
    m = size(Φ)[1]
    scale = Int(floor(m/n))
    if ( abs(Int(floor(m/n)) - m/n) > eps(1.0))
        return error
    end
    Φ2 = zeros(n)
    for i in 1:n
        Φ2[i] = 0
            for j in 1:scale
                Φ2[i] += Φ[scale*(i-1) + j] / scale
            end
    end
    return Φ2
end
function avgx(Φ, n)
    m = size(Φ)[1]
    scale = Int(floor(m/n))
    if ( abs(Int(floor(m/n)) - m/n) > eps(1.0))
        return error
    end
    nx, ny, nz = size(Φ)
    Φ2 = zeros(n, ny, nz)
    for i in 1:n
        Φ2[i,:,:] .= 0
            for j in 1:scale
                Φ2[i,:,:] .+= Φ[scale*(i-1) + j,:,:] / scale
            end
    end
    return Φ2
end
function avgy(Φ, n)
    m = size(Φ)[2]
    scale = Int(floor(m/n))
    if ( abs(Int(floor(m/n)) - m/n) > eps(1.0))
        return error
    end
    nx, ny, nz = size(Φ)
    Φ2 = zeros(nx, n, nz)
    for i in 1:n
        Φ2[:,i,:] .= 0
            for j in 1:scale
                Φ2[:,i,:] .+= Φ[:,scale*(i-1) + j,:] / scale
            end
    end
    return Φ2
end
avgxy(Φ, n) = avgx(avgy(Φ, n),n )