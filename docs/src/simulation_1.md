# Mesoscale Channel Setup
The equations of motion are [here](@ref sec:eom). The Oceananigans script is the hybrid.jl setup [here](https://github.com/sandreza/Mesoscale/blob/master/oceananigans_scripts/hybrid.jl).

# Boundary Conditions / Source Terms

- Buoyancy : Linear relaxation to the surface, insulating boundary conditions elsewhere. A sponge layer relaxation to the northern wall.
- Momentum : Linear drag for the horizontal velocities, gaussian wind-stress for the zonal velocity at the surface, stress-free boundary conditions for the meridional velocity, no penetration for the wall normal velocities.

# Table of model parameters and functional form

The default viscocity and diffusivity values for the simulation are

|   Parameter             | Value       | Units | Description |
|   :-------:             | :---:       | :---:  |:---:       |
| ``\kappa^h``           | ``0.5 \times 10^{-5} `` |  ``\frac{\text{m}^2}{\text{s}}``      | horizontal diffusivity |
| ``\kappa^v``           | ``0.5 \times 10^{-5}`` | ``\frac{\text{m}^2}{\text{s}}``           | vertical diffusivity |
| ``\nu^h``           | ``12 `` |  ``\frac{\text{m}^2}{\text{s}}``      | horizontal viscocity |
| ``\nu^v``           | ``3 \times 10^{-4}`` |   ``\frac{\text{m}^2}{\text{s}}``         | vertical viscocity |

The physical parameter values are

|   Parameter             | Value       | Units | Description |
|   :-------:             | :---:       | :---:  |:---:       |
| ``\rho^0``           | ``1024`` |  ``\frac{\text{kg}}{\text{m}^3}``      | reference density|
| ``\alpha``           | ``2 \times 10^{-4}`` | ``\frac{1}{\text{K}}``           | thermal expansion coefficient |
| ``g``           | ``9.8061 `` |  ``\frac{\text{m}}{\text{s}^2}``      | gravitational acceleration |


The domain is ``[0, L_x) \times [0, L_y] \times [-L_z, 0]``, (periodic, wall-bounded, wall-bounded) with values

|   Parameter             | Value       | Units | Description |
|   :-------:             | :---:       | :---:  |:---:       |
| ``L_x``           | ``10^3`` |  ``\text{km}``      | Zonal Length|
| ``L_y``           | ``10^{3}`` | ``\text{km}``           | Meriodional Length |
| ``L_z``           | ``3`` |  ``\text{km}``      | Depth |

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
\mathcal{S}^b(b,y,z) = - \frac{1}{\lambda^t}\left[ b - \Delta b  \frac{\exp \left( z/h \right) - \exp \left( -L_z/h \right)}{\exp \left( 0 \right) - \exp \left( -L_z/h \right)} \right] \exp \left( \frac{y - L_y}{\lambda^N} \right)
```

and the surface relaxation for buoyancy is

```math
\left. \Phi^*_b \right|_{\text{surface}} = \lambda^s \left(b - \Delta b \frac{y}{L_y} \right)
```

with parameter values

|   Parameter             | Value       | Units | Description |
|   :-------:             | :---:       | :---:  |:---:       |
| ``\lambda^s``           | ``10^{-4} `` |  ``\frac{m}{\text{s}}``      | surface buoyancy velocity constant|
| ``\Delta b``           | ``10 \times \alpha \times 10`` | ``\frac{\text{m}}{\text{ s}^2}``           | buoyancy jump|
| ``h``           | ``1`` | ``\text{km}``           | northern wall stratification e-folding length|
| ``\lambda^t``           | ``28 \times 86400`` | ``\text{s}``           | relaxation time|
| ``\lambda^N``           | ``20`` | ``\text{km}``           | northern wall horizontal e-folding length|

The surface flux for the zonal velocity, (``u`` where ``\bm{u} = (u,v,w) ``) is
```math
\left. \Phi^*_u \right|_{\text{surface}} = -\frac{\tau}{\rho_0} \left[ \exp \left( - \lambda^u \frac{(y-L_y/2)^2}{L_y^2} \right) - \exp \left( - \lambda^u \frac{(L_y/2)^2}{L_y^2} \right) \right]
```
with parameters

|   Parameter             | Value       | Units | Description |
|   :-------:             | :---:       | :---:  |:---:       |
| ``\tau``           | ``0.2 `` |  ``\frac{N}{\text{m}^2}``      | wind-stress magnitude|
| ``\lambda^u``           | ``32`` | dimensionless           | decay constant|

The linear drag was chosen to be

|   Parameter             | Value       | Units | Description |
|   :-------:             | :---:       | :---:  |:---:       |
| ``\mu``           | ``10^{-3}`` |  ``\frac{m}{\text{s}}``      | linear drag velocity parameter|

And finally the numerical parameters were

|   Parameter             | Value       | Units | Description |
|   :-------:             | :---:       | :---:  |:---:       |
| ``\Delta t``           | ``300`` |  ``\text{s}``      | time-step size|
| ``N_x``           | ``192`` |  dimensionless     | zonal grid points|
| ``N_y``           | ``192`` |  dimensionless    | meridional grid points|
| ``N_z``           | ``32`` |  dimensionless    | vertical levels|

# Comments on Parameter Choices

- ``\Delta b`` is the buoyancy jump from the bottom of the domain and the top of the domain. It is chosen so that the buoyancy should remain between ``0`` and ``\Delta b``. The relaxation on the surface is chosen so that the minimum on the surface matches with the minimum on the bottom, and similarly for the maximum. Since the smallest temperature on the surface occurs in the south, this induces isopycnals from the southern surface to the northern bottom.
- The parameters ``\lambda^N`` and ``\lambda^u`` were chosen so that there was little overlap in momentum windstress and buoyancy sponge relaxation 
- The simulation was insensitive to ``\kappa^h``, suggesting that there is a substantial amount of horizontal numerical diffusion

# Initial Condition

Firstly consider an initial condition for buoyancy of the form
```math
\begin{aligned}
b(x, y, z, t= 0 ) &= b_s(y) + b_n(z) - b_n(0)
\end{aligned} 
```
with ``b_s(L_y) = b_n(0)`` where ``b_s(y)`` is the surface relaxation profile and ``b_n(z)`` is the northern wall relaxation. Observe that this choice of initial condition implies
```math
\begin{aligned}
b(x, y = L_y, z, t= 0 ) &= b^n(z)  \\
b(x, y , z = 0, t= 0 ) &= b^s(y)  \\
\end{aligned} 
```
In the context of the simulation above this means 
```math
\begin{aligned}
b(x, y, z, t= 0 ) &=  \underbrace{\Delta b\frac{y}{L_y}}_{\text{surface relaxation}} + \underbrace{\Delta b \frac{\exp(z/h)- \exp(-L_z / h)}{1- \exp(-L_z / h)}}_{\text{northern wall relaxation}} - \underbrace{\Delta b}_{b_s(L_y) = b_n(0)} 
\end{aligned} 
```

From the initial condition for buoyancy we will derive the initial condition for velocity assuming [Geostrophic Balance](https://en.wikipedia.org/wiki/Geostrophic_current),
```math
\begin{aligned}
- (f+\beta y) v &= -\frac{1}{\rho^0}\partial_x p  \\
(f+\beta y)u &= - \frac{1}{\rho^0}\partial_y p \\ 
b &=  \frac{1}{\rho^0}\partial_z p 
\end{aligned}
```
Integrating the last equation yields
```math
p(x,y,z) - p(x,y,-L_z) = \rho^0 \int_{-L_z}^z b(y,\zeta) d\zeta,
```
which can then be differentiated with respect to ``x`` to give an expression for ``v``
```math
\begin{aligned}
(f+\beta y)v &= -\frac{1}{\rho^0}\partial_x p(x,y,-Lz) .
\end{aligned}
```
This implies that ``v`` must be independent of the vertical coordinate ``z``. In particular we assign the value ``v = 0`` at the bottom of the domain since we have a linear drag bottom boundary condition (relaxation to zero). Thus ``\partial_x p(x,y,-Lz) = 0 \Rightarrow p(x,y,-Lz) = \rho^0 F(y)`` where ``F(y)`` is some arbitrary function of ``y``.

 We now differentiate pressure with respect to ``y`` to get an expression for the zonal velocity 
 ```math
\begin{aligned}
(f+\beta y)u &= - \partial_y \left(\int_{-L_z}^z b(y,\zeta) d\zeta + F(y) \right)
 \\
&= -  \int_{-L_z}^z \partial_y b(y,\zeta) d\zeta - F'(y) 
\end{aligned}
```
Since ``\partial_y b(y,\zeta)`` is independent of ``z`` we can carry out the integral and get
 ```math
\begin{aligned}
(f+\beta y)u &= - (z+L_z)b'_s(y) - F'(y) 
\end{aligned}
```
We can determine ``F'(y) = 0`` by evaluating the above expression at ``z = -L_z`` and using that ``u  = 0`` at the bottom of the domain. We can then determine that 
 ```math
\begin{aligned}
u &= - \frac{(z+L_z)}{f+\beta y}b'_s(y)
\end{aligned}
```
The last velocity field ``w`` can be determined from the continuity equation and the condition that ``w=0`` at the surface and the bottom to get that the initial condition for ``w = 0``.

Thus in total we have the following concrete initial conditions
```math
\begin{aligned}
b^0 &= \Delta b \left(\frac{y}{L_y} +  \frac{\exp(z/h)- \exp(-L_z / h)}{1- \exp(-L_z / h)} - 1\right) \\ 
u^0 &= - \frac{(z+L_z)}{f+\beta y}\frac{\Delta b}{L_y} \\ 
v^0 &= w^0 = 0
\end{aligned}
```



