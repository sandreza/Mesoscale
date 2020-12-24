# Boiler Plate "using" in seperate file
include(pwd() * "/scripts/dependencies.jl")

# Architecture
CUDA.allowscalar(true)
arch = CPU()
FT   = Float64

# File IO 
write_slices = true
write_zonal  = true
slice_output_interval = 48hour
zonal_output_interval = 365day

time_avg_window =  zonal_output_interval / 1.0 # needs to be a float
checkpoint_interval = 365 * 1 *  day

end_time = 30day # 200 * 365day

descriptor = " "
filename = "Channel_" * descriptor

# Domain
const Lx = 1000.0kilometer 
const Ly = 1000.0kilometer
const Lz = 3.0kilometer

# Discretization
advection = WENO5()
timestepping = :RungeKutta3

# Rough resolution
Δx = Δy = 5kilometer # 5km
Δz = 100meter
# Multiple of 16 gridpoints
const Nx = round(Int, Lx / Δx / 16) # * 16
const Ny = round(Int, Ly / Δy / 16) # * 16
const Nz = round(Int, Lz / Δz / 16) * 16

# Create Grid
topology = (Periodic, Bounded, Bounded)
grid = RegularCartesianGrid(topology=topology, size=(Nx, Ny, Nz), x=(0, Lx), y=(0, Ly), z=(-Lz, 0))

# Momentum Boundary Conditions
@inline windshape(y, p) = sin( π * y / p.Ly)
@inline windstress(x, y, t, p) = - p.τ / p.ρ * windshape(y, p)
@inline ulineardrag(i, j, grid, clock, state, p) = @inbounds - p.μ * state.u[i, j, 1]
@inline vlineardrag(i, j, grid, clock, state, p) = @inbounds - p.μ * state.v[i, j, 1]
# Zonal Velocity
top_u_bc = BoundaryCondition(Flux, windstress, parameters = bc_params)
bottom_u_bc =  BoundaryCondition(Flux, ulineardrag, discrete_form = true, parameters = bc_params)
u_bcs = UVelocityBoundaryConditions(grid, top = top_u_bc, bottom = bottom_u_bc)
# Meridional Velocity
bottom_v_bc =  BoundaryCondition(Flux, vlineardrag, discrete_form = true, parameters = bc_params)
v_bcs = VVelocityBoundaryConditions(grid, bottom = bottom_v_bc)

# Buoyancy Boundary Conditions Forcing. Note: Flux convention opposite of Abernathy
@inline cutoff(j, grid, p ) = grid.yC[j] > p.Qᵇ_cutoff ? -0.0 : 1.0
@inline surface_flux(j, grid, p) = p.Qᵇ * cos(3π * grid.yC[j] / p.Ly) * cutoff(j, grid, p)
@inline relaxation(i, j, grid, clock, state, p) = @inbounds surface_flux(j, grid, p)
top_b_bc = BoundaryCondition(Flux, relaxation, discrete_form = true, parameters = bc_params)
b_bcs = TracerBoundaryConditions(grid, top = top_b_bc)


# Sponge layers
relu(y) = (abs(y) + y) * 0.5
# Northern Wall Relaxation
@inline relaxation_profile_north(k, grid, p) = p.ΔB * ( exp(grid.zC[k]/p.h) - exp(-p.Lz/p.h) ) / (1 - exp(-p.Lz/p.h))
function Fb_function(i, j, k, grid, clock, state, p)
    return @inbounds - (1/p.λᵗ)  * (state.b[i,j,k] 
        - relaxation_profile_north(k, grid, p)) * relu( (grid.yC[j]-p.Lsponge) / (p.Ly - p.Lsponge))
end
Fb = Forcing(Fb_function, parameters = bc_params, discrete_form = true)
forcings = (b = Fb, ) 