using Oceananigans
using Oceananigans.Grids
using JLD2
import Oceananigans.Fields: interpolate
using Oceananigans.Fields: fractional_indices, _interpolate

"""
    interpolate(field, LX, LY, LZ, grid, x, y, z)
Interpolate `field` to the physical point `(x, y, z)` using trilinear interpolation. The location of
the field is specified with `(LX, LY, LZ)` and the field is defined on `grid`.
Note that this is a lower-level `interpolate` method defined for use in CPU/GPU kernels.
"""
@inline function interpolate(field, LX, LY, LZ, grid, x, y, z)
    i, j, k = fractional_indices(x, y, z, (LX, LY, LZ), grid)

    # We use mod and trunc as CUDA.modf is not defined.
    # For why we use Base.unsafe_trunc instead of trunc see:
    # https://github.com/CliMA/Oceananigans.jl/issues/828
    # https://github.com/CliMA/Oceananigans.jl/pull/997
    ξ, i = mod(i, 1), Base.unsafe_trunc(Int, i)
    η, j = mod(j, 1), Base.unsafe_trunc(Int, j)
    ζ, k = mod(k, 1), Base.unsafe_trunc(Int, k)

    return _interpolate(field, ξ, η, ζ, i+1, j+1, k+1)
end

function interpolate(field::AbstractArray, LX, LY, LZ, grid, x::AbstractArray, y::AbstractArray, z::AbstractArray)
    newgrid = zeros((length(x), length(y), length(z)))
    for k in eachindex(z), j in eachindex(y), i in eachindex(x)
        newgrid[i,j,k] = interpolate(field, LX, LY, LZ, grid, x[i], y[j], z[k])
    end
    return newgrid
end

function getgrid(file)
    coarsetobe = jldopen(file)
    grid = coarsetobe["grid"]
    close(coarsetobe)
    return grid
end

function interpolate(states::Array{Array{S, 3},1}, grid, newgrid) where {S}
    newx, newy, newz = nodes((Cell, Cell, Cell), new_grid, reshape=true)
    interpolatedstates = [interpolate(states, LX, LY, LZ, grid, newx, newy, newz) for state in states]
    return interpolatedstates
end
##
#=
# Example
files = [pwd() * "/Channel_1_checkpoint_iteration404905.jld2",
         pwd() * "/Channel_2_checkpoint_iteration728850.jld2",
         pwd() * "/Channel_3_checkpoint_iteration1087211.jld2",
         pwd() * "/Channel_4_checkpoint_iteration6160221.jld2",
         pwd() * "/Channel_8_checkpoint_iteration3164301.jld2",
         pwd() * "/Channel_16_checkpoint_iteration6317902.jld2",
         pwd() * "/Channel_24_checkpoint_iteration8612852.jld2", 
         pwd() * "/Channel_32_checkpoint_iteration1272125.jld2"
]

file = files[end]
coarsetobe = jldopen(file)
grid = coarsetobe["grid"]
u = coarsetobe["velocities"]["u"]
LX, LY, LZ = typeof.(u["location"])
array = u["data"]
# newpoint = interpolate(array, LX, LY, LZ, grid, 0.0, 0.0, 00.0)
i, j, k = fractional_indices(0, 0, -2000, (LX, LY, LZ), grid)
##

topo = (Periodic, Bounded, Bounded)
domain = (x=(0, 1e6), y=(0, 1e6), z=(-3000, 0))
new_grid = RegularCartesianGrid(topology=topo, size=(150, 150, 12); domain...)
newx, newy, newz = nodes((Cell, Cell, Cell), new_grid, reshape=true)

field = u["data"]
newgrid = interpolate(field, LX, LY, LZ, grid, newx, newy, newz)

# newgridvals = interpolate.(Ref(array), Ref(LX), Ref(LY), Ref(LZ), Ref(grid),newx, newy, newx)

scene = volume(0..1, 0..1, 0..1, newgrid)
volume!(scene, 1..2, 0..1, 0..1, field)
=#
