using Printf, Statistics
using Oceananigans
using Oceananigans.BoundaryConditions
using Oceananigans.Fields
using Oceananigans.OutputWriters
using Oceananigans.Diagnostics
using Oceananigans.Utils
using Oceananigans.AbstractOperations
using CUDA
CUDA.allowscalar(true)
arch = GPU()
FT   = Float64

write_output = false
geostrophic_balance = false
output_interval = 365 * 24hour # 48hour makes nice movies
time_avg_window =  output_interval / 2.0 # needs to be a float
checkpoint_interval = 365 * 1 * day

end_time = 100 * 365day
const scale = 20;
filename_1 = "Ridge_" * string(scale)

const Lx = scale * 50kilometer # 1000km 
const Ly = scale * 50kilometer # 1000km
const Lz = 3kilometer          # 3km

# Rough resolution
Δx = Δy = scale * 250meter # 5km
Δz = 100meter

const Nx = round(Int, Lx/ Δx / 16) * 16
const Ny = round(Int, Ly/ Δy / 16) * 16
const Nz = round(Int, Lz/ Δz / 16) * 16

topology = (Periodic, Bounded, Bounded)
grid = RegularCartesianGrid(topology=topology, size=(Nx, Ny, Nz), x=(0, Lx), y=(0, Ly), z=(-Lz, 0))

const f = -1e-4
const β = 1 * 10^(-11)
coriolis = FPlane(FT, f=f)
coriolis = BetaPlane(FT, f₀ = f, β = β)

α = 2e-4  # Thermal expansion coefficient [K⁻¹]
eos = LinearEquationOfState(FT, α=α, β=0)
buoyancy = BuoyancyTracer()

κh = 0.5e-5
νh = 12.0 
κv = 0.5e-5 
νv = 3e-4 
closure = AnisotropicDiffusivity(νx = νh, νy = νh, νz =νv, κx = κh, κy = κh, κz=κv)

bc_params = (
    Lx = Lx,
    Ly = Ly,
    τ = 0.2, # [N m⁻²] Zonal stress
    ρ = 1024, # [kg / m³]
    μ = 1.0e-3,  # [m/s] linear drag decay scale
    H = Lz, 
    λˢ = 1e-4, # [m/s] relaxation boundary condition
    h = 1000.0, # [m] relexaction profile scale
    ΔB = 10 * 2e-3, # buoyancy jump
    Lz = Lz,
    λᵘ = 32, # surface forcing e-folding length scale
    λᵗ = 7.0 * 86400.0, # [s]
    λᴺ = 2.0 * 10^4, #[m] northern wall e-folding scale
    τᵇ = 600.0, # [s] bottom sponge relaxation time-scale (basically 2Δt)
    λᵇ = 200.0 # [m] bottom sponge relaxation length
)

const h = 1000.0 # [m]
const ΔB = 10 * 2e-3
const ΔT = 10 * 2e-3


@inline wind_shape(y, p) = exp( -(y - p.Ly/2)^2 / (p.Ly^2 / p.λᵘ) ) - exp( -( p.Ly/2)^2 / (p.Ly^2 / p.λᵘ) )
@inline wind_stress(x, y, t, p) = - p.τ / p.ρ * wind_shape(y, p)

@inline relaxation_profile(j, grid, p) = p.ΔB * (grid.yC[j]/ p.Ly)
@inline relaxation(i, j, grid, clock, state, p) = @inbounds p.λˢ * ( state.tracers.b[i, j, grid.Nz] - relaxation_profile(j, grid, p))

## Sponge layers
# U, V bottom relaxation
const ridge_height = 1800.0 #[m]
const ridge_base = -2800.0 #[m]
@inline ridge_shape(x,z,L) = -z + ridge_base + ridge_height * exp(-40 *(x-L/2)^2 / L^2)
@inline smoothed_ridge(x,z,L) =  (tanh(ridge_shape(x,z,L)) +1)/2
@inline ridge(x,z,L) = O < ridge_shape(x,z,L) ? 1.0 : -0.0
@inline function Fu_function(i, j, k, grid, clock, state, p)
    return @inbounds ( -1.0/p.τᵇ * state.velocities.u[i,j,k] *
                    smoothed_ridge(grid.xF[i], grid.zC[k], p.Lx) 
                      )
end

@inline function Fv_function(i, j, k, grid, clock, state, p)
    return @inbounds ( -1.0/p.τᵇ * state.velocities.v[i,j,k] *
                       smoothed_ridge(grid.xC[i], grid.zC[k], p.Lx)
                      )
end

@inline function Fw_function(i, j, k, grid, clock, state, p)
    return @inbounds ( -1.0/p.τᵇ * state.velocities.w[i,j,k] *
                       smoothed_ridge(grid.xC[i], grid.zF[k], p.Lx)
                      )
end

# Northern Wall Relaxation
@inline relaxation_profile_north(k, grid, p) = p.ΔB * ( exp(grid.zC[k]/p.h) - exp(-p.Lz/p.h) ) / (1 - exp(-p.Lz/p.h))
function Fb_function(i, j, k, grid, clock, state, p)
    return @inbounds - (1/p.λᵗ)  * (state.tracers.b[i,j,k] 
        - relaxation_profile_north(k, grid, p)) * exp( (grid.yC[j] - p.Ly) / p.λᴺ ) 
end

Fu = ParameterizedForcing(Fu_function, bc_params)
Fv = ParameterizedForcing(Fv_function, bc_params)
Fw = ParameterizedForcing(Fw_function, bc_params)
Fb = ParameterizedForcing(Fb_function, bc_params)
forcings = ModelForcing(u = Fu, v = Fv, w = Fw, b = Fb)

# Boundary Conditions
# Buoyancy
top_b_bc = BoundaryCondition(Flux, relaxation, discrete_form = true, parameters = bc_params)
b_bcs = TracerBoundaryConditions(grid, top = top_b_bc)
# Zonal Velocity
top_u_bc = BoundaryCondition(Flux, wind_stress, parameters = bc_params)
u_bcs = UVelocityBoundaryConditions(grid, top = top_u_bc)

# Initial Conditions
ε(σ) = σ * randn()
if !geostrophic_balance
    # initial condition the same as the northern wall relaxation
    B₀(x, y, z) = ΔB * ( exp(z/h) - exp(-Lz/h) ) / (1 - exp(-Lz/h)) + ε(1e-8)
else
    U₀(x,y,z) = (z + Lz)/(f + β * y) * (ΔB / Ly)
    B₀(x, y, z) = ΔB * (y / Ly +  ( exp(z/h) - exp(-Lz/h) ) / (1 - exp(-Lz/h)) - 1)
end

# boundary conditions
bcs = (b = b_bcs,  u = u_bcs)

## checkpointing / model construction
include(pwd() * "/oceananigans_scripts/checkpointing.jl")

## Diagnostics
include(pwd() * "/oceananigans_scripts/diagnostics.jl")

## Set timestep
# Δt is defined in checkpointing
Δt_wizard = TimeStepWizard(cfl = 0.3, Δt = Δt, max_change = 1.1, max_Δt = 300.0)
cfl = AdvectiveCFL(Δt_wizard)

## Progress Printing
include(pwd() * "/oceananigans_scripts/progress_printer.jl")

simulation = Simulation(model, Δt=Δt_wizard, stop_time=end_time, progress=print_progress, iteration_interval=Ni)

include(pwd() * "/oceananigans_scripts/output_writing.jl")

## Run Simulation
run!(simulation)
