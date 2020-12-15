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
arch = CPU()
FT   = Float64

write_output = true 
geostrophic_balance = false
output_interval = 365 * 24hour # 48hour makes nice movies
time_avg_window = floor(Int, output_interval / 2)
checkpoint_interval = 365 * 4 * day

end_time = 100*365day
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
    λᵗ = 28.0 * 86400.0, # [s]
    λᴺ = 2.0 * 10^4, #[m] northern wall e-folding scale
)

# Initial Conditions
U₀(x,y,z) = (z + Lz)/(f + β * y) * (ΔB / Ly)
B₀(x, y, z) = ΔB * (y / Ly +  ( exp(z/h) - exp(-Lz/h) ) / (1 - exp(-Lz/h)) - 1)

## checkpointing
model = IncompressibleModel(
            architecture = arch,
                float_type = FT,
                    grid = grid,
                coriolis = coriolis,
                buoyancy = buoyancy,
                closure = closure,
                tracers = (:b,),
)
set!(model, u = U₀, b=B₀)

## Visualization
using NetCDF, Plots, GLMakie, AbstractPlotting
using Printf, Statistics
using ImageTransformations, Colors, JLD2
using AbstractPlotting.MakieLayout

statenames = ("u", "v", "w", "b")
function visualize(model::Oceananigans.AbstractModel)
    # fix me
    u = Array(interior(model.velocities.u))
    v = Array(interior(model.velocities.v))
    w = Array(interior(model.velocities.w))
    b = Array(interior(model.tracers.b))
    # fix above
    states = [u, v, w, b]
    statenames = ("u", "v", "w", "b")
    visualize(states, statenames = statenames)
    return nothing
end

function visualize(states; statenames = string.(1:length(states)))
    stateindex = collect(1:length(states))
    statenode = Node(stateindex[4])
    colorchoices = [:balance, :thermal, :dense, :deep, :curl, :thermometer]
    colornode = Node(colorchoices[1])

    state = @lift(states[$statenode])
    scene, layout = layoutscene()
    lscene = layout[1:4, 2:4] = LScene(scene)

    @lift begin
        x = size(states[$statenode])[1] 
        y = size(states[$statenode])[2] 
        z = size(states[$statenode])[3] 
    end

    clims = @lift(extrema(states[$statenode])) 

    cmap_rgb = @lift(to_colormap($colornode))
    titlename = @lift(" "^10 * " Field " * statenames[$statenode] * " "^10) # use padding and appropriate centering
    supertitle = layout[1,2] = LText(scene, titlename , textsize = 50, color = :black)

    volume!(lscene, 0..x, 0..y, 0..z, state, 
            camera = cam3d!, 
            colormap = cmap_rgb, 
            colorrange = clims)

    statemenu = LMenu(scene, options = zip(statenames, stateindex))
    colormenu = LMenu(scene, options = zip(colorchoices, colorchoices))

    on(statemenu.selection) do s
        statenode[] = s
    end

    on(colormenu.selection) do s
        colornode[] = s
    end

    layout[1, 1] = vgrid!(
        LText(scene, "State", width = nothing),
        statemenu,
        LText(scene, "Color", width = nothing),
        colormenu,
    )
    display(scene)
    return nothing
end

