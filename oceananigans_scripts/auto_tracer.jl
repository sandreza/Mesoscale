using Printf
using Statistics

using Oceananigans
using Oceananigans.Units
using Oceananigans.OutputReaders: FieldTimeSeries
using Oceananigans.Grids: xnode, ynode, znode
using Random
using JLD2

import Oceananigans.Grids: ynode

hydrostatic = false
stretched_grid = false
load_from_checkpoint = true
const Lx = 1000kilometers # zonal domain length [m]
const Ly = 2000kilometers # meridional domain length [m]


symbol_list = Symbol[]
case = "trial"
if case == "trial"
    jlist = [1,2]
    klist = [1,2,3,4]
end
for j in jlist, k in klist
    push!(symbol_list, Meta.parse("c_k"*string(j) * "_j"*string(k)))
end
t_list = Tuple(symbol_list)

# file output
prefix = "relaxation_channel_tracers_restarted_"* "smooth_forcing" * "_case_" * case

# time
Δt = 300 * 2
stop_time = 200 * 365days # 300years
# number of grid points
Nx = 16 * 8
Ny = Nx * 2
Nz = 32

# stretched grid 
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

# Tracer Forcings
function_name_list = Symbol[]
jk_list = []
for j in jlist, k in klist
    push!(function_name_list, Meta.parse("forcing_c_k"*string(j) * "_j"*string(k)))
    push!(jk_list, (j,k))
end
f_n_list = Tuple(function_name_list)
f_list = []
for i in eachindex(f_n_list)
    f_n = f_n_list[i]
    jj, kk = jk_list[i]
    @eval begin
        @inline function $f_n(i, j, k, grid, clock, model_fields, parameters) 
            y = ynode(Center(), j, grid)
            z = znode(Center(), k, grid)
            forcing = cos( acos(π * (-z / parameters.Lz) ) * $kk) * cos(acos(π * (y / parameters.Ly) )* $jj)
            return forcing
        end
        push!($f_list, Forcing($f_n, discrete_form = true, parameters = parameters))
    end
end

tracer_name_list = (t_list..., :b)
tracer_function_list = (Tuple(f_list)..., Fb)
tracer_forcing_list = (; zip(tracer_name_list, tracer_function_list)...)


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

model = NonhydrostaticModel(architecture = arch,
                                    grid = grid,
                                    advection = WENO5(),
                                    buoyancy = BuoyancyTracer(),
                                    coriolis = coriolis,
                                    closure = (closure),
                                    tracers = tracer_name_list,
                                    boundary_conditions = (b=b_bcs, u=u_bcs, v=v_bcs),
                                    forcing = tracer_forcing_list,
                                    )


@info "Built $model."

#####
##### Initial conditions
#####

# resting initial condition
Random.seed!(1)
ε(σ) = σ * randn()
bᵢ(x, y, z) = parameters.ΔB * ( exp(z/parameters.h) - exp(-Lz/parameters.h) ) / (1 - exp(-Lz/parameters.h)) + ε(1e-8)

# set model initial conditions
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

if load_from_checkpoint 
    jl_file = jldopen("relaxation_channel_nh_iteration15768000.jld2", "r+")
    u = Oceananigans.CUDA.CuArray(jl_file["velocities"]["u"]["data"])
    v = Oceananigans.CUDA.CuArray(jl_file["velocities"]["v"]["data"])
    w = Oceananigans.CUDA.CuArray(jl_file["velocities"]["w"]["data"])

    b = Oceananigans.CUDA.CuArray(jl_file["tracers"]["b"]["data"])
    #=
    c1 = Oceananigans.CUDA.CuArray(jl_file["tracers"]["c1"]["data"])
    c2 = Oceananigans.CUDA.CuArray(jl_file["tracers"]["c2"]["data"])
    c3 = Oceananigans.CUDA.CuArray(jl_file["tracers"]["c3"]["data"])
    c4 = Oceananigans.CUDA.CuArray(jl_file["tracers"]["c4"]["data"])
    =#

    model.velocities.u.data[:] .= u[:]
    model.velocities.v.data[:] .= v[:]
    model.velocities.w.data[:] .= w[:]
    model.tracers.b.data[:]    .= b[:]

else
    nothing
end

simulation = Simulation(model, Δt=wizard, stop_time=stop_time, progress=print_progress, iteration_interval=1000)

#####
##### Diagnostics
#####

u, v, w = model.velocities
b = model.tracers.b

# ζ = ComputedField(∂x(v) - ∂y(u))

averaged_outputs = Dict(
    :u  => AveragedField(u,  dims=(1,)),
    :v  => AveragedField(v,  dims=(1,)),
    :w  => AveragedField(w,  dims=(1,)),

    :uu => AveragedField(u * u, dims= 1),
    :vv => AveragedField(v * v, dims= 1),
    :ww => AveragedField(w * w, dims= 1),
    :uv => AveragedField(u * v, dims= 1),
    :vw => AveragedField(w * v, dims= 1),
    :uw => AveragedField(w * u, dims= 1),
)

for i in eachindex(tracer_name_list)
    push!(averaged_outputs, tracer_name_list[i] => AveragedField(getproperty(model.tracers, tracer_name_list[i]), dims=(1,))) 
    push!(averaged_outputs, tracer_name_list[i] => AveragedField(v * getproperty(model.tracers, tracer_name_list[i]), dims=(1,))) 
    push!(averaged_outputs, tracer_name_list[i] => AveragedField(w * getproperty(model.tracers, tracer_name_list[i]), dims=(1,))) 
end

simulation.output_writers[:averages] = JLD2OutputWriter(model, averaged_outputs,
                                                    schedule = AveragedTimeInterval(10*365days, window=10*365days, stride=10),
                                                    prefix = prefix * "_averages",
                                                    verbose = true,
                                                    force = true)

 @info "Running the simulation..."

#  try
tic = time()
run!(simulation, pickup=false)
toc = time()