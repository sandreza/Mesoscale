# Mesoscale Channel Setup
The equations of motion are [here](@ref sec:eom). The Oceananigans script is the hybrid.jl setup [here](https://github.com/sandreza/Mesoscale/blob/master/oceananigans_scripts/hybrid.jl).

# Boundary Conditions / Source Terms

- Buoyancy : Flux boundary conditions on the surface, insulating boundary conditions elsewhere. A sponge layer relaxation to the northern wall.
- Momentum : Linear drag for the horizontal velocities, sinusoidal wind-stress for the zonal velocity at the surface, stress-free boundary conditions for the meridional velocity, no penetration for the wall normal velocities.

# Numerics 

- Advection Scheme: 3rd-5th order, WENO5
- Timestepping: 1st-2nd order, Quasi-AB2

# Table of model parameters and functional form

The default viscocity and diffusivity values for the simulation are

|   Parameter             | Value       | Units | Description |
|   :-------:             | :---:       | :---:  |:---:       |
| ``\kappa^h``           | ``0.0 `` |  ``\frac{\text{m}^2}{\text{s}}``      | horizontal diffusivity |
| ``\kappa^v``           | ``0.5 \times 10^{-5}`` | ``\frac{\text{m}^2}{\text{s}}``           | vertical diffusivity |
| ``\nu^h``           | ``12 `` |  ``\frac{\text{m}^2}{\text{s}}``      | horizontal viscocity |
| ``\nu^v``           | ``3 \times 10^{-4}`` |   ``\frac{\text{m}^2}{\text{s}}``         | vertical viscocity |

The physical parameter values are

|   Parameter             | Value       | Units | Description |
|   :-------:             | :---:       | :---:  |:---:       |
| ``\rho^0``           | ``999.8`` |  ``\frac{\text{kg}}{\text{m}^3}``      | reference density|
| ``\alpha``           | ``2 \times 10^{-4}`` | ``\frac{1}{\text{K}}``           | thermal expansion coefficient |
| ``g``           | ``9.8061 `` |  ``\frac{\text{m}}{\text{s}^2}``      | gravitational acceleration |
| ``c^p``           | ``3994 `` |  ``\frac{\text{J}}{\text{K} \text{kg}}``      | heat capacity |


The domain is ``[0, L_x) \times [0, L_y] \times [-L_z, 0]``, (periodic, wall-bounded, wall-bounded) with values

|   Parameter             | Value       | Units | Description |
|   :-------:             | :---:       | :---:  |:---:       |
| ``L_x``           | ``10^3`` |  ``\text{km}``      | Zonal Length|
| ``L_y``           | ``2 \times 10^{3}`` | ``\text{km}``           | Meriodional Length |
| ``L_z``           | ``2985`` |  ``\text{m}``      | Depth |

We use a ``\beta``-plane coriolis 
```math
\bm{f} = (f + \beta y)\hat{z}
```

with ``f`` and ``\beta`` values

|   Parameter             | Value       | Units | Description |
|   :-------:             | :---:       | :---:  |:---:       |
| ``f``           | ``-10^{-4} `` |  ``\frac{1}{\text{s}}``      | constant term|
| ``\beta``           | ``10^{-11}`` | ``\frac{1}{\text{m}} \frac{1}{\text{ s}}``           | linear term|

The sponge layer for buoyancy at northern wall is of the form

```math
\mathcal{S}^b(b,y,z) = - \frac{1}{\lambda^t}\left[ b - \Delta b  \frac{\exp \left( z/h \right) - \exp \left( -L_z/h \right)}{\exp \left( 0 \right) - \exp \left( -L_z/h \right)} \right] \text{relu}\left( \frac{y - L_{\text{sponge}}}{L_y - L_{\text{sponge}}} \right)
```
where the [relu](https://en.wikipedia.org/wiki/Rectifier_(neural_networks)) function is the rectified linear unit,

and the surface forcing for buoyancy is

```math
\left. \Phi^*_b \right|_{\text{surface}} = 
\begin{cases}
\frac{Q \alpha g}{\rho^0 c^p} \cos \left( \frac{3 \pi}{L_y} y \right) &  \text{ for } y \leq \frac{5}{6} L_y \\ 
0 & \text{otherwise}
\end{cases}
```

with parameter values

|   Parameter             | Value       | Units | Description |
|   :-------:             | :---:       | :---:  |:---:       |
| ``Q``           | ``10 `` |  ``\frac{W}{\text{m}^2}``      | surface flux|
| ``h``           | ``1`` | ``\text{km}``           | northern wall stratification e-folding length|
| ``\lambda^t``           | ``7 \times 86400`` | ``\text{s}``           | relaxation time|
| ``L_{\text{sponge}}``           | ``1980`` | ``\text{km}``           | northern wall meridional activation level |

The surface flux for the zonal velocity, (``u`` where ``\bm{u} = (u,v,w) ``) is
```math
\left. \Phi^*_u \right|_{\text{surface}} = -\frac{\tau}{\rho_0} \sin \left( \frac{\pi y}{L_y} \right)
```
with parameters

|   Parameter             | Value       | Units | Description |
|   :-------:             | :---:       | :---:  |:---:       |
| ``\tau``           | ``0.2 `` |  ``\frac{N}{\text{m}^2}``      | wind-stress magnitude|

The linear drag was chosen to be

|   Parameter             | Value       | Units | Description |
|   :-------:             | :---:       | :---:  |:---:       |
| ``\mu``           | ``10^{-3}`` |  ``\frac{m}{\text{s}}``      | linear drag velocity parameter|

And finally the numerical parameters were

|   Parameter             | Value       | Units | Description |
|   :-------:             | :---:       | :---:  |:---:       |
| ``\Delta t``           | ``300`` |  ``\text{s}``      | time-step size|
| ``N_x``           | ``384`` |  dimensionless     | zonal grid points|
| ``N_y``           | ``192`` |  dimensionless    | meridional grid points|
| ``N_z``           | ``32`` |  dimensionless    | vertical levels|

# Comments on Parameter Choices

- The simulation was insensitive to ``\kappa^h``, suggesting that there is a substantial amount of horizontal numerical diffusion, potentially in every direction

# Initial Condition
The velocity field was chosen to be zero and the buoyancy field starts as the northern wall relaxation profile (meaning no lateral gradients). A small amount of noise was added to the buoyancy field to aid in the transition to turublence