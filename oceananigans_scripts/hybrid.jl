using Printf
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
checkpoint_interval = 365 * 5 * day

end_time = 100 * 365day
const scale = 20;
filename_1 = "Hybrid_" * string(scale)

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
)

const h = 1000.0 # [m]
const ΔB = 10 * 2e-3

@inline wind_shape(y, p) = exp( -(y - p.Ly/2)^2 / (p.Ly^2 / p.λᵘ) ) - exp( -( p.Ly/2)^2 / (p.Ly^2 / p.λᵘ) )
@inline wind_stress(x, y, t, p) = - p.τ / p.ρ * wind_shape(y, p)
@inline τ₁₃_linear_drag(i, j, grid, clock, state, p) = @inbounds - p.μ * state.velocities.u[i, j, 1]
@inline τ₂₃_linear_drag(i, j, grid, clock, state, p) = @inbounds - p.μ * state.velocities.v[i, j, 1]



@inline relaxation_profile(j, grid, p) = p.ΔB * (grid.yC[j]/ p.Ly)
@inline relaxation(i, j, grid, clock, state, p) = @inbounds p.λˢ * ( state.tracers.b[i, j, grid.Nz] - relaxation_profile(j, grid, p))

# Sponge layers
# Northern Wall Relaxation
@inline relaxation_profile_north(k, grid, p) = p.ΔB * ( exp(grid.zC[k]/p.h) - exp(-p.Lz/p.h) ) / (1 - exp(-p.Lz/p.h))
function Fb_function(i, j, k, grid, clock, state, p)
    return @inbounds - (1/p.λᵗ)  * (state.tracers.b[i,j,k] 
        - relaxation_profile_north(k, grid, p)) * exp( (grid.yC[j] - p.Ly) / p.λᴺ ) 
end

Fb = ParameterizedForcing(Fb_function, bc_params)
forcings = ModelForcing(b = Fb)

# Boundary Conditions
# Buoyancy
top_b_bc = BoundaryCondition(Flux, relaxation, discrete_form = true, parameters = bc_params)
b_bcs = TracerBoundaryConditions(grid, top = top_b_bc)
# Zonal Velocity
top_u_bc = BoundaryCondition(Flux, wind_stress, parameters = bc_params)
bottom_u_bc =  BoundaryCondition(Flux, τ₁₃_linear_drag, discrete_form = true, parameters = bc_params)
u_bcs = UVelocityBoundaryConditions(grid, top = top_u_bc, bottom = bottom_u_bc)
# Meridional Velocity
bottom_v_bc =  BoundaryCondition(Flux, τ₂₃_linear_drag, discrete_form = true, parameters = bc_params)
v_bcs = VVelocityBoundaryConditions(grid, bottom = bottom_v_bc)

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
bcs = (b = b_bcs,  u = u_bcs, v = v_bcs)

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