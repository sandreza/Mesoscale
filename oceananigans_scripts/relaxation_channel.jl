using Printf, Revise
using Statistics

using Oceananigans
using Oceananigans.Units
using Oceananigans.OutputReaders: FieldTimeSeries
using Oceananigans.Grids: xnode, ynode, znode
using Random

import Oceananigans.Grids: ynode

hydrostatic = false
stretched_grid = false
const Lx = 1000kilometers # zonal domain length [m]
const Ly = 2000kilometers # meridional domain length [m]

# file output
prefix = "relaxation_channel_nh"

# time
Δt = 300 * 2
stop_time = 300years
# number of grid points
Nx = 16 * 8
Ny = Nx * 2
Nz = 32

# stretched grid 
#=
k_center = collect(1:Nz)
Δz_center = @. 10 * 1.104^(Nz - k_center)
const Lz = sum(Δz_center)
z_faces = vcat([-Lz], -Lz .+ cumsum(Δz_center))
z_faces[Nz+1] = 0
=#
const Lz = 3000.00

println("Largest Lz is ", Lz)

arch = GPU()
FT = Float64
topology = (Periodic, Bounded, Bounded)
if stretched_grid
    grid = VerticallyStretchedRectilinearGrid(architecture = arch,
                                            topology = topology,
                                            size = (Nx, Ny, Nz),
                                            halo = (3, 3, 3),
                                            x = (0, Lx),
                                            y = (0, Ly),
                                            z_faces = z_faces)
else
    grid = RegularRectilinearGrid(topology=topology, 
                            size=(Nx, Ny, Nz), 
                            x=(0, Lx), y=(0, Ly), z=(-Lz, 0))
end



@info "Built a grid: $grid."



#####
##### Boundary conditions
#####

α  = 2e-4     # [K⁻¹] thermal expansion coefficient 
g  = 9.8061   # [m/s²] gravitational constant
cᵖ = 3994.0   # [J/K]  heat capacity
ρ  = 999.8    # [kg/m³] reference density


parameters = (
    Ly = Ly,  
    Lz = Lz,    
    Qᵇ = 10/(ρ * cᵖ) * α * g,            # buoyancy flux magnitude [m² s⁻³]    
    y_shutoff = 5/6 * Ly,                # shutoff location for buoyancy flux [m]
    τ = 0.2,                            # surface kinematic wind stress [kg m⁻¹ s⁻²]
    μ = 1.1e-3,                          # linear drag [m/s]  
    ΔB = 8 * α * g,                      # surface vertical buoyancy gradient [s⁻²]
    H = Lz,                              # domain depth [m]
    h = 1000.0,                          # exponential decay scale of stable stratification [m]
    y_sponge = 19/20 * Ly,               # southern boundary of sponge layer [m]
    λt = 7.0days,                         # relaxation time scale [s]
    α  = 2e-4,     # [K⁻¹] thermal expansion coefficient 
    g  = 9.8061,   # [m/s²] gravitational constant
    cᵖ = 3994.0,   # [J/K]  heat capacity
    ρ  = 999.8,   # [kg/m³] reference density
    λᵘ = 32.0, # surface forcing e-folding length scale
    λˢ = 1e-4, # [m/s] relaxation boundary condition
)


@inline relaxation_profile(y, p) = p.ΔB * (y/ p.Ly)
@inline relaxation(i, j, grid, clock, state, p) = @inbounds p.λˢ * ( state.b[i, j, grid.Nz] - relaxation_profile(y, p))
# ifelse(y < p.y_shutoff, p.Qᵇ * cos(3π * y / p.Ly), 0.0)
@inline function buoyancy_flux(i, j, grid, clock, state, p)
    y = ynode(Center(), j, grid)
    return @inbounds p.λˢ * ( state.b[i, j, grid.Nz] - relaxation_profile(y, p))
end

buoyancy_flux_bc = FluxBoundaryCondition(buoyancy_flux, discrete_form=true, parameters=parameters)

@inline wind_shape(y, p) = exp( -(y - p.Ly/2)^2 / (p.Ly^2 / p.λᵘ) ) - exp( -( p.Ly/2)^2 / (p.Ly^2 / p.λᵘ) ) 
@inline wind_stress(y, p) = - p.τ / p.ρ * wind_shape(y, p)
# - p.τ / p.ρ * sin(π * y / p.Ly)
@inline function u_stress(i, j, grid, clock, model_fields, p)
    y = ynode(Center(), j, grid)
    return wind_stress(y, p)
end

u_stress_bc = FluxBoundaryCondition(u_stress, discrete_form=true, parameters=parameters)



@inline u_drag(i, j, grid, clock, model_fields, p) = @inbounds - p.μ * model_fields.u[i, j, 1] 
@inline v_drag(i, j, grid, clock, model_fields, p) = @inbounds - p.μ * model_fields.v[i, j, 1]

u_drag_bc = FluxBoundaryCondition(u_drag, discrete_form=true, parameters=parameters)
v_drag_bc = FluxBoundaryCondition(v_drag, discrete_form=true, parameters=parameters)

b_bcs = FieldBoundaryConditions(top = buoyancy_flux_bc)

u_bcs = FieldBoundaryConditions(top = u_stress_bc, bottom = u_drag_bc)
v_bcs = FieldBoundaryConditions(bottom = v_drag_bc)


#####
##### Coriolis
#####

const f = -1e-4
const β = 1 * 10^(-11)
coriolis = BetaPlane(FT, f₀ = f, β = β)

#####
##### Forcing and initial condition
#####

@inline initial_buoyancy(z, p) = p.ΔB * (exp(z / p.h) - exp(-p.Lz / p.h)) / (1 - exp(-p.Lz / p.h))
@inline mask(y, p) = max(0.0, y - p.y_sponge) / (p.Ly - p.y_sponge)


@inline function buoyancy_relaxation(i, j, k, grid, clock, model_fields, p)
    timescale = p.λt
    y = ynode(Center(), j, grid)
    z = znode(Center(), k, grid)
    target_b = initial_buoyancy(z, p)
    b = @inbounds model_fields.b[i, j, k]
    return - 1 / timescale  * mask(y, p) * (b - target_b)
end

Fb = Forcing(buoyancy_relaxation, discrete_form = true, parameters = parameters)

# closure

κh = 0.5e-5 # [m²/s] horizontal diffusivity
νh = 12.0   # [m²/s] horizontal viscocity
κz = 0.5e-5 # [m²/s] vertical diffusivity
νz = 3e-4   # [m²/s] vertical viscocity

closure = AnisotropicDiffusivity(νh=νh, νz=νz, κh=κh, κz=κz)


convective_adjustment = ConvectiveAdjustmentVerticalDiffusivity(convective_κz = 1.0,
                                                               convective_νz = 0.0)


#####
##### Model building
#####

@info "Building a model..."

if hydrostatic
    model = HydrostaticFreeSurfaceModel(architecture = arch,
                                        grid = grid,
                                        free_surface = ImplicitFreeSurface(),
                                        momentum_advection = WENO5(),
                                        tracer_advection = WENO5(),
                                        buoyancy = BuoyancyTracer(),
                                        coriolis = coriolis,
                                        closure = (closure, convective_adjustment),
                                        tracers = :b,
                                        boundary_conditions = (b=b_bcs, u=u_bcs, v=v_bcs),
                                        forcing = (b=Fb,),
                                        )
else
    model = NonhydrostaticModel(architecture = arch,
                                        grid = grid,
                                        advection = WENO5(),
                                        buoyancy = BuoyancyTracer(),
                                        coriolis = coriolis,
                                        closure = (closure),
                                        tracers = :b,
                                        boundary_conditions = (b=b_bcs, u=u_bcs, v=v_bcs),
                                        forcing = (b=Fb,),
                                        )
end

@info "Built $model."

#####
##### Initial conditions
#####

# resting initial condition
Random.seed!(1)
ε(σ) = σ * randn()
bᵢ(x, y, z) = parameters.ΔB * ( exp(z/parameters.h) - exp(-Lz/parameters.h) ) / (1 - exp(-Lz/parameters.h)) + ε(1e-8)

set!(model, b=bᵢ)

#####
##### Simulation building

wizard = Δt # TimeStepWizard(cfl=0.1, Δt=Δt, max_change=1.0, max_Δt=Δt)

wall_clock = [time_ns()]

function print_progress(sim)
    @printf("[%05.2f%%] i: %d, t: %s, wall time: %s, max(u): (%6.3e, %6.3e, %6.3e) m/s, next Δt: %s\n",
            100 * (sim.model.clock.time / sim.stop_time),
            sim.model.clock.iteration,
            prettytime(sim.model.clock.time),
            prettytime(1e-9 * (time_ns() - wall_clock[1])),
            maximum(abs, sim.model.velocities.u),
            maximum(abs, sim.model.velocities.v),
            maximum(abs, sim.model.velocities.w),
 #           prettytime(sim.Δt.Δt))
            prettytime(sim.Δt))

    wall_clock[1] = time_ns()
    
    return nothing
end

simulation = Simulation(model, Δt=wizard, stop_time=stop_time, progress=print_progress, iteration_interval=1000)

#####
##### Diagnostics
#####

u, v, w = model.velocities
b = model.tracers.b

# ζ = ComputedField(∂x(v) - ∂y(u))

averaged_outputs = Dict(
    :u => AveragedField(u, dims=(1,)),
    :v => AveragedField(v, dims=(1,)),
    :w => AveragedField(w, dims=(1,)),
    :b => AveragedField(b, dims=(1,)),

    :uu => AveragedField(u * u, dims= 1),
    :vv => AveragedField(v * v, dims= 1),
    :ww => AveragedField(w * w, dims= 1),
    :uv => AveragedField(u * v, dims= 1),
    :vw => AveragedField(w * v, dims= 1),
    :uw => AveragedField(w * u, dims= 1),
    :ub => AveragedField(u * b, dims= 1),
    :vb => AveragedField(v * b, dims= 1),
    :wb => AveragedField(w * b, dims= 1),
)

# #####
# ##### Build checkpointer and output writer
# #####

simulation.output_writers[:checkpointer] = Checkpointer(model,
                                                    schedule = TimeInterval(5 * 365days),
                                                    prefix = prefix,
                                                    force = true)

simulation.output_writers[:averages] = JLD2OutputWriter(model, averaged_outputs,
                                                    schedule = AveragedTimeInterval(10*365days, window=10*365days, stride=10),
                                                    prefix = prefix * "_averages",
                                                    verbose = true,
                                                    force = true)

 @info "Running the simulation..."

#  try
run!(simulation, pickup=false)