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

write_output = true 
output_interval = 60 * 24hour # 48hour
const scale = 20;
filename_1 = "old_setup_" * string(scale)

const Lx = scale * 50kilometer # 2000km cos( 5 * π * y / p.Ly ) 
const Ly = scale * 50kilometer # 2000km
const Lz = 3kilometer          # 3km

Δx = Δy = scale * 250meter # 8x larger?
Δz = 100meter

const Nx = round(Int, Lx/ Δx / 16) * 16
const Ny = round(Int, Ly/ Δy / 16) * 16
const Nz = round(Int, Lz/ Δz / 16) * 16

topology = (Periodic, Bounded, Bounded)
grid = RegularCartesianGrid(topology=topology, size=(Nx, Ny, Nz), x=(0, Lx), y=(0, Ly), z=(-Lz, 0))

f = -1e-4
β = 1 * 10^(-11)
coriolis = FPlane(FT, f=f)
coriolis = BetaPlane(FT, f₀ = f, β = β)

α = 2e-4  # Thermal expansion coefficient [K⁻¹]
eos = LinearEquationOfState(FT, α=α, β=0)
buoyancy = BuoyancyTracer()
# buoyancy = SeawaterBuoyancy(FT, equation_of_state=eos, constant_salinity=true)
const vis_scale = 1.0;
κh = 0.5e-5*vis_scale
νh = 12.0   # Horizontal diffusivity and viscosity [m²/s]
κv = 0.5e-5 *vis_scale
νv = 3e-4 * vis_scale  # Vertical diffusivity and viscosity [m²/s]
closure = AnisotropicDiffusivity(νx = νh, νy = νh, νz =νv, κx = κh, κy = κh, κz=κv)
# closure = AnisotropicMinimumDissipation(FT)

bc_params = (
    Ly = Ly,
    B½ = 1.96e-7 / 40,    # Buoyancy flux at midchannel [m²/s³]
    Lᶠ = scale * 2kilometer, # Characteristic length scale of the forcing [m]
    τ = 0.2, # [N m⁻²] Zonal stress
    ρ = 1024, # [kg / m³]
    μ = 1.0e-3,  # [m/s] linear drag decay scale
    H = Lz, 
    λ = 1e-4, # [m/s] relaxation boundary condition
    h = 1000, # [m] relexaction profile scale
    ΔB = 10 * 2e-3, # buoyancy jump
    Lz = Lz
)

@inline wind_stress(x, y, t, p) = - p.τ / p.ρ * sin( π*y / p.Ly)
@inline buoyancy_flux(x, y, t, p) = -p.B½ * y / p.Ly

@inline τ₁₃_linear_drag(i, j, grid, clock, state, p) = @inbounds -p.μ * state.velocities.u[i, j, 1]
@inline τ₂₃_linear_drag(i, j, grid, clock, state, p) = @inbounds -p.μ * state.velocities.v[i, j, 1]

const h = 1000.0 # [m]
const ΔT = 10 * 2e-3
@inline function zC(k)
    return - Lz + (k-0.5) / Nz * Lz
end

@inline function yC(j)
    return  (j-0.5) / Nx * Ly
end
@inline relaxation_profile(j, p) = p.ΔB * (yC(j)/ p.Ly)
@inline relaxation(i, j, grid, clock, state, p) = @inbounds p.λ * ( state.tracers.b[i, j, grid.Nz] - relaxation_profile(j, p))

@inline relaxation_profile_north(k, p) = ΔT * ( exp(zC(k)/h) - exp(-Lz/h) ) / (1 - exp(-Lz/h))
@inline relaxation_north(i, k, grid, clock, state, p) = @inbounds 1e3 * p.λ * ( state.tracers.b[i, grid.Ny, k] - relaxation_profile_north(k, p))

# Buoyancy
top_b_bc = ParameterizedBoundaryCondition(Flux, relaxation, bc_params)
north_b_bc = ParameterizedBoundaryCondition(Flux, relaxation_north, bc_params)
b_bcs = TracerBoundaryConditions(grid, top = top_b_bc, north = north_b_bc)

# Meridonal Velocity
bottom_v_bc =  ParameterizedBoundaryCondition(Flux, τ₂₃_linear_drag, bc_params)
v_bcs = VVelocityBoundaryConditions(grid, bottom = bottom_v_bc)

# Zonal Velocity
u_velocity_flux_bf = BoundaryFunction{:z, Face, Cell}(wind_stress, bc_params)
top_u_bc = FluxBoundaryCondition(u_velocity_flux_bf)
bottom_u_bc =  ParameterizedBoundaryCondition(Flux, τ₁₃_linear_drag, bc_params)
u_bcs = UVelocityBoundaryConditions(grid, top = top_u_bc, bottom = bottom_u_bc)

ε(σ) = σ * randn()
# make initial condition same as relaxation to northern wall
B₀(x, y, z) = ΔT * ( exp(z/h) - exp(-Lz/h) ) / (1 - exp(-Lz/h)) + ε(1e-8)
# boundary conditions
bcs = (b = b_bcs,  u = u_bcs, v = v_bcs)

# checkpointing
searchdir(path, key) = filter(x -> occursin(key, x), readdir(path))
checkpoints = searchdir(pwd(), filename_1 * "_checkpoint_iteration")
if length(checkpoints) > 0
    checkpointed = true
    checkpoint_filepath = joinpath(pwd(), checkpoints[end])
    @info "Restoring from checkpoint: $checkpoint_filepath"
    model = restore_from_checkpoint(checkpoint_filepath, boundary_conditions = bcs)
else
    checkpointed = false
	model = IncompressibleModel(
	           architecture = arch,
	             float_type = FT,
	                   grid = grid,
	               coriolis = coriolis,
	               buoyancy = buoyancy,
	                closure = closure,
	                tracers = (:b,),
	    boundary_conditions = bcs
	)
end
if !checkpointed
	set!(model, b=B₀)
end

## End Model description

fields = Dict(
    "u" => model.velocities.u,
    "v" => model.velocities.v,
    "w" => model.velocities.w,
    "b" => model.tracers.b
)

surface_output_writer =
    NetCDFOutputWriter(model, fields, filename= filename_1 * "_surface.nc",
			           time_interval=output_interval, zC=Nz, zF=Nz)

middepth_output_writer =
    NetCDFOutputWriter(model, fields, filename= filename_1 * "_middepth.nc",
                       time_interval=output_interval, zC=Int(Nz/2), zF=Int(Nz/2))

zonal_output_writer =
    NetCDFOutputWriter(model, fields, filename= filename_1 * "_zonal.nc",
                       time_interval=output_interval, yC=Int(Ny/2), yF=Int(Ny/2))

meridional_output_writer =
    NetCDFOutputWriter(model, fields, filename= filename_1 * "_meridional.nc",
                       time_interval=output_interval, xC=Int(Nx/2), xF=Int(Nx/2))

checkpointer = Checkpointer(model, prefix = filename_1 * "_checkpoint", time_interval = 30day, force = true)
###
# kinetic_energy = @at (Cell, Cell, Cell) (u^2 + v^2 + w^2) / 2 where u = models.velocities.u, etc
#bouyancy profile
Uz = Average(model.velocities.u; return_type=Array, dims = (1,))
Bz = Average(model.tracers.b; return_type=Array, dims = (1,))
# Example = ZonalAverage(∂x(model.tracers.b) * ∂y(model.tracers.u); return_type = Array)
# Create output writer that writes vertical profiles to JLD2 output files.
zonal_averages = Dict(
	"Uz" => model -> Uz(model)[1, 2:end-1, 2:end-1],
	"Bz" => model -> Bz(model)[1, 2:end-1, 2:end-1],
)

output_attributes = Dict(
    "Uz" => Dict("longname" => "Zonal Average Velocity in the x-direction", "units" => "m/s"),
    "Bz" => Dict("longname" => "Zonal Average Buoyancy", "units" => "m/s²"),
)
# Should probably output error if this is not supplied for nonfield objects
dimensions = Dict(
	"Uz" => ("yC", "zC"),
	"Bz" => ("yC", "zC"),
)

zonal_average_output_writer = NetCDFOutputWriter(model, zonal_averages, filename =  filename_1 * "_zonal_average.nc", time_interval=output_interval, output_attributes=output_attributes, dimensions = dimensions)
# close(zonal_average_output_writer)
###
if !checkpointed
	Δt= 120.0
else
	Δt = 300.0
end
Δt_wizard = TimeStepWizard(cfl=0.3, Δt = Δt, max_change=1.1, max_Δt= 300.0)
cfl = AdvectiveCFL(Δt_wizard)


# Take Ni "intermediate" time steps at a time before printing a progress
# statement and updating the time step.
Ni = 1000

function print_progress(simulation)
    model = simulation.model
    i, t = model.clock.iteration, model.clock.time

    progress = 100 * (model.clock.time / end_time)

    umax = maximum(abs, model.velocities.u.data.parent)
    vmax = maximum(abs, model.velocities.v.data.parent)
    wmax = maximum(abs, model.velocities.w.data.parent)

    @printf("[%05.2f%%] i: %d, t: %.2e days, umax: (%6.3e, %6.3e, %6.3e) m/s, CFL: %6.4e, next Δt: %.1e s\n",
    	    progress, i, t / day, umax, vmax, wmax, cfl(model), Δt_wizard.Δt)
end

end_time = 103*365day
simulation = Simulation(model, Δt=Δt_wizard, stop_time=end_time, progress=print_progress, iteration_interval=Ni)
if write_output
    simulation.output_writers[:surface] = surface_output_writer
    simulation.output_writers[:middepth] = middepth_output_writer
    simulation.output_writers[:zonal] = zonal_output_writer
    simulation.output_writers[:meridional] = meridional_output_writer
    simulation.output_writers[:zonal_average] = zonal_average_output_writer
end
simulation.output_writers[:checkpoint] = checkpointer
###
run!(simulation)
write_output(model, checkpointer)


