using Oceananigans
using Oceananigans.Models: ShallowWaterModel

Lx, Ly, Lz = 2π, 2π, 1
Nx, Ny = 64*4*4, 64*4*4

grid = RectilinearGrid(size = (Nx, Ny),
                       x = (0, Lx), y = (-Ly/2, Ly/2),
                       topology = (Periodic, Periodic, Flat))

const U = 1.0         # Maximum jet velocity

f = 1           # Coriolis parameter
g = 10       # Gravitational acceleration

model = ShallowWaterModel(architecture = GPU(),
                          timestepper = :RungeKutta3,
                          advection = WENO5(),
                          grid = grid,
                          gravitational_acceleration = g,
                          )
h̄(x, y, z) = Lz
ū(x, y, z) = U * sech(y)^2
ω̄(x, y, z) = 2 * U * sech(y)^2 * tanh(y)

small_amplitude = 1e-4

 uⁱ(x, y, z) = ū(x, y, z) + small_amplitude * exp(-y^2) * randn()
uhⁱ(x, y, z) = uⁱ(x, y, z) * h̄(x, y, z)

ū̄h(x, y, z) = ū(x, y, z) * h̄(x, y, z)

set!(model, uh = ū̄h, h = h̄)

uh, vh, h = model.solution

# Build velocities
u = uh / h
v = vh / h

# Build and compute mean vorticity discretely
ω = ComputedField(∂x(v) - ∂y(u))
compute!(ω)

# Copy mean vorticity to a new field
ωⁱ = Field(Face, Face, Nothing, model.architecture, model.grid)
ωⁱ .= ω

# Use this new field to compute the perturbation vorticity
ω′ = ComputedField(ω - ωⁱ)

set!(model, uh = uhⁱ)

simulation = Simulation(model, Δt = 5e-4, stop_time = 200)
tic = time()
run!(simulation)
toc = time()
println("the time for the simulation is ", toc - tic)