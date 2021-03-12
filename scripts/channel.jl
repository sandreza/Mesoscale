# Boiler Plate "using" in seperate file
include(pwd() * "/scripts/dependencies.jl")

# Architecture
CUDA.allowscalar(true)
arch = CPU()
FT   = Float64

# File IO 
ic_load = true

write_slices = false
write_zonal  = false
slice_output_interval = 48hour
zonal_output_interval = 365day

time_avg_window =  zonal_output_interval / 2  # needs to be a float
checkpoint_interval = 365 * 10 *  day

resolution = 1
descriptor = string(resolution)
filename = "Channel_" * descriptor

if ic_load
    filepath = pwd() * "/" * getlatest(filename)
end
## Domain
const Lx = 1000.0kilometer 
const Ly = 1000.0kilometer
const Lz = 3.0kilometer

# Discretization
maxΔt = 1100.0 * 16 / resolution  # [s]
Δt =  ic_load ? maxΔt : Δt = 300.0 * 16 / resolution # [s]

end_time = 200 * 365day 
advection   = WENO5()
timestepper = :RungeKutta3
# Rough target resolution
Δx = Δy = 5kilometer # 5km
Δz = 100meter
# Multiple of 16 gridpoints
const Nx = round(Int, Lx / Δx / 16) * resolution
const Ny = round(Int, Ly / Δy / 16) * resolution
const Nz = round(Int, Lz / Δz / 16) * 16
# Create Grid
topology = (Periodic, Bounded, Bounded)
grid = RegularCartesianGrid(topology=topology, 
                            size=(Nx, Ny, Nz), 
                            x=(0, Lx), y=(0, Ly), z=(-Lz, 0))

# Parameters
const f = -1e-4
const β = 1 * 10^(-11)
coriolis = FPlane(FT, f=f)
coriolis = BetaPlane(FT, f₀ = f, β = β)

α  = 2e-4     # [K⁻¹] Thermal expansion coefficient 
g  = 9.8061   # [m/s²] gravitational constant
cᵖ = 3994.0   # [J/K]  heat capacity
ρ  = 999.8    # [kg/m³] density
const h = 1000.0     # [m] e-folding length scale for northern sponge
const ΔB = 8 * α * g # [m/s²] total change in buoyancy from surface to bottom
eos = LinearEquationOfState(FT, α=α, β=0)
buoyancy = BuoyancyTracer()

κh = 0.5e-5 # [m²/s] horizontal diffusivity
νh = 12.0   # [m²/s] horizontal viscocity
κv = 0.5e-5 # [m²/s] vertical diffusivity
νv = 3e-4   # [m²/s] vertical viscocity

closure = AnisotropicDiffusivity(νx = νh, νy = νh, νz =νv, 
                                 κx = κh, κy = κh, κz=κv)

parameters = (
    Ly = Ly,                   # y-domain length
    τ = 0.2,                   # [N m⁻²] Zonal stress
    ρ = ρ,                     # [kg / m³]
    μ = 1.1e-3,                # [m/s]  linear drag
    H = Lz,                    # [m]
    h = h,                     # [m]    relexaction profile scale
    ΔB = ΔB,                   # [m/s²] buoyancy jump
    Lz = Lz,                   # [m]
    Lsponge = 900kilometer,   # [m]
    λᵗ = 7.0day,               # [s]
    Qᵇ = 10/(ρ * cᵖ) * α * g,  # [m² / s³]
    Qᵇ_cutoff = Ly * 5/6.      # [m]
)

# Momentum Boundary Conditions
@inline windshape(y, p) = sin( π * y / p.Ly)
@inline windstress(x, y, t, p) = - p.τ / p.ρ * windshape(y, p)
@inline ulineardrag(i, j, grid, clock, state, p) = @inbounds - p.μ * state.u[i, j, 1]
@inline vlineardrag(i, j, grid, clock, state, p) = @inbounds - p.μ * state.v[i, j, 1]
# Zonal Velocity
top_u_bc = BoundaryCondition(Flux, windstress, parameters = parameters)
bottom_u_bc =  BoundaryCondition(Flux, ulineardrag, discrete_form = true, parameters = parameters)
u_bcs = UVelocityBoundaryConditions(grid, top = top_u_bc, bottom = bottom_u_bc)
# Meridional Velocity
bottom_v_bc =  BoundaryCondition(Flux, vlineardrag, discrete_form = true, parameters = parameters)
v_bcs = VVelocityBoundaryConditions(grid, bottom = bottom_v_bc)

# Buoyancy Boundary Conditions Forcing. Note: Flux convention opposite of Abernathy
@inline cutoff(j, grid, p ) = grid.yC[j] > p.Qᵇ_cutoff ? -0.0 : 1.0
@inline surface_flux(j, grid, p) = p.Qᵇ * cos(3π * grid.yC[j] / p.Ly) * cutoff(j, grid, p)
@inline relaxation(i, j, grid, clock, state, p) = @inbounds surface_flux(j, grid, p)
top_b_bc = BoundaryCondition(Flux, relaxation, discrete_form = true, parameters = parameters)
b_bcs = TracerBoundaryConditions(grid, top = top_b_bc)

# Save boundary conditions as named tuple
bcs = (b = b_bcs,  u = u_bcs, v = v_bcs,)

# Forcing Functions
# Sponge layers
relu(y) = (abs(y) + y) * 0.5
# Northern Wall Relaxation
@inline relaxation_profile_north(k, grid, p) = p.ΔB * ( exp(grid.zC[k]/p.h) - exp(-p.Lz/p.h) ) / (1 - exp(-p.Lz/p.h))
function Fb_function(i, j, k, grid, clock, state, p)
    return @inbounds - (1/p.λᵗ)  * (state.b[i,j,k] 
        - relaxation_profile_north(k, grid, p)) * relu( (grid.yC[j]-p.Lsponge) / (p.Ly - p.Lsponge))
end
Fb = Forcing(Fb_function, parameters = parameters, discrete_form = true)

# Record forcings
forcings = (b = Fb, ) 

# Model Setup
model = IncompressibleModel(
           architecture = arch,
             float_type = FT,
                   grid = grid,
               coriolis = coriolis,
               buoyancy = buoyancy,
                closure = closure,
                tracers = (:b,),
    boundary_conditions = bcs,
                forcing = forcings,
              advection = advection,
            timestepper = timestepper,
)

# Initial Conditions
if ic_load
    println("loading from " * filepath)
    ic!(model, filepath, ArrayType = archarray(arch))
else
    ε(σ) = σ * randn()
    B₀(x, y, z) = ΔB * ( exp(z/h) - exp(-Lz/h) ) / (1 - exp(-Lz/h)) + ε(1e-8)
    set!(model, b=B₀)
end

# Input Output Details
include(pwd() * "/scripts/channel_io.jl")

# Create Simulation
Δt_wizard = TimeStepWizard(cfl = 1.0, Δt = Δt, max_change = 1.05, max_Δt = maxΔt)
cfl = AdvectiveCFL(Δt_wizard)
simulation = Simulation(model, Δt=Δt_wizard, 
                        stop_time=end_time,
                        iteration_interval=Ni, 
                        progress=print_progress)
# IO
if write_slices
    simulation.output_writers[:surface] = surface_output_writer
    simulation.output_writers[:middepth] = middepth_output_writer
    simulation.output_writers[:zonal] = zonal_output_writer
    simulation.output_writers[:meridional] = meridional_output_writer
end
if write_zonal
    simulation.output_writers[:statistics] = zonal_statistics
end
simulation.output_writers[:checkpoint] = checkpointer

## Run
run!(simulation)

