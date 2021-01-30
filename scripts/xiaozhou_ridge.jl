include(pwd() * "/scripts/dependencies.jl")

advection   = WENO5()
timestepper = :RungeKutta3
arch = CPU()
FT = Float64
const Nx = 200
const Ny = 1
const Nz = 300

const Lx = 1e6  # [m]
const Ly = 1e6
const Lz = 1500.0 # [m]
topology = (Periodic, Periodic, Bounded)
grid = RegularCartesianGrid(topology=topology, 
                            size=(Nx, Ny, Nz), 
                            x=(0, Lx), y=(0,Ly), z=(-Lz, 0))

const f = 1e-4
coriolis = FPlane(FT, f=f)

α  = 2e-4     # [K⁻¹] Thermal expansion coefficient 
g  = 9.81     # [m/s²] gravitational constant
eos = LinearEquationOfState(FT, α=α, β=0)
buoyancy = BuoyancyTracer()

κh = 1e-3 # [m²/s] horizontal diffusivity
νh = 1e-3   # [m²/s] horizontal viscocity
κv = 1e-3 # [m²/s] vertical diffusivity
νv = 1e-3   # [m²/s] vertical viscocity

closure = AnisotropicDiffusivity(νx = νh, νy = νh, νz =νv, 
                                 κx = κh, κy = κh, κz=κv)

# Immersed Boundary
const ridge_height = 400.0 #[m]
const ridge_base = -1000.0  #[m]
@inline ridge_shape(x, z, L) = -z + ridge_base + ridge_height * sin(2π / Lx * x)
@inline ridge(x,y,z) = 0 < ridge_shape(x, z, 1e6)

model = IncompressibleModel(
           architecture = arch,
             float_type = FT,
                   grid = grid,
               coriolis = coriolis,
               buoyancy = buoyancy,
                closure = closure,
                tracers = (:b,),
              advection = advection,
            timestepper = timestepper,
            immersed_boundary = ridge,
)


N² = 1e-6
U₀(x,y,z) = 1e-3 * randn()
V₀(x,y,z) = 1e-3 * randn()
B₀(x,y, z) = N² * (z + Lz)
set!(model, b=B₀, u = U₀, v = V₀)


Δt = 200.0
maxΔt = 200.0
end_time = 10day
Ni = 1000
Δt_wizard = TimeStepWizard(cfl = 1.0, Δt = Δt, max_change = 1.05, max_Δt = maxΔt)
cfl = AdvectiveCFL(Δt_wizard)

function print_progress(simulation)
    model = simulation.model
    i, t = model.clock.iteration, model.clock.time

    progress = 100 * (model.clock.time / end_time)
    @printf("i: %d, t: %.2e days ", i, t / day)
    println("-----------------------------")
end

simulation = Simulation(model, Δt=Δt_wizard, 
                        stop_time=end_time,
                        iteration_interval=Ni, 
                        progress=print_progress)
##
run!(simulation)
##
vstates = [Array(interior(model.velocities[velocity]))[:,1,:] for velocity in keys(model.velocities)]
vstatenames = [string(velocity) for velocity in keys(model.velocities)]
tstates = [Array(interior(model.tracers[tracer])[:,1,:]) for tracer in keys(model.tracers)]
tstatenames = [string(tracer) for tracer in keys(model.tracers)]
states = vcat(vstates, tstates)
statenames = vcat(vstatenames, tstatenames)
scene = visualize(states, statenames = statenames)