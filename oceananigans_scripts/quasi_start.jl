using Printf
using Oceananigans
using Oceananigans.BoundaryConditions
using Oceananigans.Fields
using Oceananigans.OutputWriters
using Oceananigans.Diagnostics
using Oceananigans.Utils
using Oceananigans.AbstractOperations

arch = GPU()
FT   = Float64

write_output = true 
geostrophic_balance = true
output_interval = 2day # 2day makes nice movies
checkpointer_output_time = 5day
end_time =  30day
const scale = 20;
filename_1 = "Geostrophic_Start_" * string(scale)

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
    λᵗ = 28.0 * 86400.0, # [s]
    λᴺ = 2.0 * 10^4, #[m] northern wall e-folding scale
)

@inline wind_stress(x, y, t, p) = - p.τ / p.ρ * ( exp( -(y - p.Ly/2)^2 / (p.Ly^2 / p.λᵘ) ) - exp( -( p.Ly/2)^2 / (p.Ly^2 / p.λᵘ) ) )
@inline τ₁₃_linear_drag(i, j, grid, clock, state, p) = @inbounds - p.μ * state.velocities.u[i, j, 1]
@inline τ₂₃_linear_drag(i, j, grid, clock, state, p) = @inbounds - p.μ * state.velocities.v[i, j, 1]

const h = 1000.0 # [m]
const ΔB = 10 * 2e-3

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
if !geostrophic_balance
    # initial condition the same as the northern wall relaxation
    B₀(x, y, z) = ΔB * ( exp(z/h) - exp(-Lz/h) ) / (1 - exp(-Lz/h)) + ε(1e-8)
else
    U₀(x,y,z) = -(z + Lz)/(f + β * y) * (ΔB / Ly)
    B₀(x, y, z) = ΔB * (y / Ly +  ( exp(z/h) - exp(-Lz/h) ) / (1 - exp(-Lz/h)) - 1)
end

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
    print("Setting Initial Conditions")
    if !geostrophic_balance
        set!(model, b=B₀)
    else
        set!(model, u = U₀, b=B₀)
    end
end

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

checkpointer = Checkpointer(model, prefix = filename_1 * "_checkpoint", time_interval = checkpointer_output_time, force = true)
##
#bouyancy profile
Uz = Average(model.velocities.u; return_type=Array, dims = (1,))
Bz = Average(model.tracers.b; return_type=Array, dims = (1,))

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

##
if !checkpointed
	Δt = 300.0
else
	Δt = 300.0
end
Δt_wizard = TimeStepWizard(cfl=0.3, Δt= Δt, max_change=1.1, max_Δt= 300.0)
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


simulation = Simulation(model, Δt=Δt_wizard, stop_time=end_time, progress=print_progress, iteration_interval=Ni)
if write_output
    simulation.output_writers[:surface] = surface_output_writer
    simulation.output_writers[:middepth] = middepth_output_writer
    simulation.output_writers[:zonal] = zonal_output_writer
    simulation.output_writers[:meridional] = meridional_output_writer
    simulation.output_writers[:zonal_average] = zonal_average_output_writer
end
simulation.output_writers[:checkpoint] = checkpointer
##
run!(simulation)