using Printf, Statistics
using Oceananigans
using Oceananigans.BoundaryConditions
using Oceananigans.Fields
using Oceananigans.OutputWriters
using Oceananigans.Diagnostics
using Oceananigans.Utils
using Oceananigans.AbstractOperations
using Oceananigans.Advection
using CUDA
CUDA.allowscalar(true)
arch = CPU()
FT   = Float64

write_slices = false
write_zonal  = true
geostrophic_balance = false
advection_scheme = WENO5()
timestepping_scheme = :RungeKutta3
slice_output_interval = 48hour
zonal_output_interval = 365day
time_avg_window =  zonal_output_interval / 1.0 # needs to be a float
checkpoint_interval = 365 * 1 *  day

end_time = 0.1day # 200 * 365day
const scale = 20;
filename_1 = "NewTrial_" * string(scale)

const Lx = scale * 50.0kilometer  # 1000km 
const Ly = scale * 50.0kilometer  # 1000km
const Lz = 2985.0                 # 3km

# Rough resolution
Δx = Δy = scale * 250meter # 5km
Δz = 100meter

const Nx = round(Int, Lx/ Δx / 16) # * 16
const Ny = round(Int, Ly/ Δy / 16) # * 16
const Nz = round(Int, Lz/ Δz / 16) * 16

topology = (Periodic, Bounded, Bounded)
grid = RegularCartesianGrid(topology=topology, size=(Nx, Ny, Nz), x=(0, Lx), y=(0, Ly), z=(-Lz, 0))

const f = -1e-4
const β = 1 * 10^(-11)
coriolis = FPlane(FT, f=f)
coriolis = BetaPlane(FT, f₀ = f, β = β)

α = 2e-4     # Thermal expansion coefficient [K⁻¹]
g = 9.8061   # gravitational constant
cᵖ = 3994.0  # heat capacity
ρ = 999.8    # density

const h = 1000.0 # [m]
const ΔB = 8 * α * g

eos = LinearEquationOfState(FT, α=α, β=0)
buoyancy = BuoyancyTracer()

κh = 0.0
νh = 12.0 
κv = 0.5e-5 
νv = 3e-4 
closure = AnisotropicDiffusivity(νx = νh, νy = νh, νz =νv, κx = κh, κy = κh, κz=κv)

bc_params = (
    Ly = Ly,
    τ = 0.2,                   # [N m⁻²] Zonal stress
    ρ = ρ,                     # [kg / m³]
    μ = 1.1e-3,                # [m/s]  linear drag
    H = Lz,                    # [m]
    h = h,                # [m]    relexaction profile scale
    ΔB = ΔB,          # [m/s²] buoyancy jump
    Lz = Lz,                   # [m]
    Lsponge = 1980kilometer,   # [m]
    λᵗ = 7.0*86400.0,          # [s]
    Qᵇ = 10/(ρ * cᵖ) * α * g,  # [m² / s³]
    Qᵇ_cutoff = Ly * 5/6.      # [m]
)

# Momentum Boundary Conditions
@inline wind_shape(y, p) = sin( π * y / p.Ly)
@inline wind_stress(x, y, t, p) = - p.τ / p.ρ * wind_shape(y, p)
@inline τ₁₃_linear_drag(i, j, grid, clock, state, p) = @inbounds - p.μ * state.u[i, j, 1]
@inline τ₂₃_linear_drag(i, j, grid, clock, state, p) = @inbounds - p.μ * state.v[i, j, 1]

# Buoyancy Boundary Conditions Forcing Note: Flux convention opposite of Abernathy
@inline cutoff(j, grid, p ) = grid.yC[j] > p.Qᵇ_cutoff ? -0.0 : 1.0
@inline surface_flux(j, grid, p) = p.Qᵇ * cos(3π * grid.yC[j] / p.Ly) * cutoff(j, grid, p)
@inline relaxation(i, j, grid, clock, state, p) = @inbounds surface_flux(j, grid, p)

# Sponge layers
relu(y) = (abs(y) + y) * 0.5
# Northern Wall Relaxation
@inline relaxation_profile_north(k, grid, p) = p.ΔB * ( exp(grid.zC[k]/p.h) - exp(-p.Lz/p.h) ) / (1 - exp(-p.Lz/p.h))
function Fb_function(i, j, k, grid, clock, state, p)
    return @inbounds - (1/p.λᵗ)  * (state.b[i,j,k] 
        - relaxation_profile_north(k, grid, p)) * relu( (grid.yC[j]-p.Lsponge) / (p.Ly - p.Lsponge))
end
Fb_function(i, j, k, grid, clock, state, p) = state.b[i,j,k] 
Fb = Forcing(Fb_function, parameters = bc_params, discrete_form = true)
forcings = (b = Fb, )

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
B₀(x, y, z) = ΔB * ( exp(z/h) - exp(-Lz/h) ) / (1 - exp(-Lz/h)) + ε(1e-8)

# boundary conditions
bcs = (b = b_bcs,  u = u_bcs, v = v_bcs,)

# model setup
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
              advection = advection_scheme,
            timestepper = timestepping_scheme,
)
set!(model, b=B₀)

## Progress printer
Ni = 100

function print_progress(simulation)
    model = simulation.model
    i, t = model.clock.iteration, model.clock.time

    progress = 100 * (model.clock.time / end_time)

    umax = maximum(abs, model.velocities.u.data.parent)
    vmax = maximum(abs, model.velocities.v.data.parent)
    wmax = maximum(abs, model.velocities.w.data.parent)
    meanb = mean(view(model.tracers.b.data, 2:Nx-1, 2:Ny-1, 2:Nz-1))

    @printf("[%05.2f%%] i: %d, t: %.2e days, umax: (%6.3e, %6.3e, %6.3e) m/s, CFL: %6.4e, next Δt: %.1e s\n",
            progress, i, t / day, umax, vmax, wmax, cfl(model), Δt_wizard.Δt)
    println(" ")
    @printf("The mean value of b is %.2e", meanb)
    println(" ")
end

## Diagnostics
fields = Dict(
    "u" => model.velocities.u,
    "v" => model.velocities.v,
    "w" => model.velocities.w,
    "b" => model.tracers.b
)

surface_output_writer =
    NetCDFOutputWriter(model, fields, filepath = filename_1 * "_surface.nc",
			           schedule=TimeInterval(slice_output_interval), field_slicer = FieldSlicer(k = Nz))

middepth_output_writer =
    NetCDFOutputWriter(model, fields, filepath = filename_1 * "_middepth.nc",
                       schedule=TimeInterval(slice_output_interval), field_slicer = FieldSlicer(k = Int(floor(Nz/2))))

zonal_output_writer =
    NetCDFOutputWriter(model, fields, filepath = filename_1 * "_zonal.nc",
                       schedule=TimeInterval(slice_output_interval), field_slicer = FieldSlicer(j = Int(floor(Ny/2))))

meridional_output_writer =
    NetCDFOutputWriter(model, fields, filepath = filename_1 * "_meridional.nc",
                       schedule=TimeInterval(slice_output_interval), field_slicer = FieldSlicer(i = Int(floor(Nx/2))))
function debug_f(debug)
    if debug
        close(surface_output_writer)
        close(middepth_output_writer)
        close(zonal_output_writer)
        close(meridional_output_writer)
    end
    return nothing
end


## Zonal and Time Averages
u, v, w = model.velocities
b = model.tracers.b

u_scratch = XFaceField(model.architecture, model.grid)
v_scratch = YFaceField(model.architecture, model.grid)
w_scratch = ZFaceField(model.architecture, model.grid)
b_scratch =  CellField(model.architecture, model.grid)

zonal_fields = Dict(
    :u => AveragedField( u, dims=(1,)),
    :v => AveragedField( v, dims=(1,)),
    :w => AveragedField( w, dims=(1,)),
    :b => AveragedField( b, dims=(1,)),
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

zonal_statistics = JLD2OutputWriter(model, zonal_fields,
                     time_averaging_window = time_avg_window,
                             time_interval = zonal_output_interval,
                                    prefix = filename_1 * "_zonal_averages",
                                     force = true)

## Checkpointer
checkpointer = Checkpointer(model, prefix = filename_1 * "_checkpoint", time_interval = checkpoint_interval, force = true)

## Create the Simulation
Δt = 300.0
Δt_wizard = TimeStepWizard(cfl = 1.0, Δt = Δt, max_change = 1.05, max_Δt = 300.0)
cfl = AdvectiveCFL(Δt_wizard)
simulation = Simulation(model, Δt=Δt_wizard, stop_time=end_time, iteration_interval=Ni, progress=print_progress)
## Output Writing
simulation.output_writers[:surface] = surface_output_writer
simulation.output_writers[:middepth] = middepth_output_writer
simulation.output_writers[:zonal] = zonal_output_writer
simulation.output_writers[:meridional] = meridional_output_writer
simulation.output_writers[:statistics] = zonal_statistics
simulation.output_writers[:checkpoint] = checkpointer
## Run the Model
run!(simulation)