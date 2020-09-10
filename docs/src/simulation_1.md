# Mesoscale Channel Setup
The equations of motion are [here](@ref sec:eom). 

# Boundary Conditions / Source Terms

- Buoyancy : Linear relaxation to the surface, insulating boundary conditions elsewhere. A sponge layer relaxation to the northern wall.
- Momentum : Linear drag and no penetration on the bottom, gaussian wind-stress for the zonal velocity at the surface, stress-free boundary conditions for the meridional velocity and no penetration at the surface

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
| ``\alpha``           | ``2 \times 10^{-4}`` | ``\frac{\text{m}}{\text{s}^2 \text{K}}``           | thermal expansion coefficient |
| ``g``           | ``9.8061 `` |  ``\frac{\text{m}}{\text{s}^2}``      | gravitational acceleration |


The domain is ``[0, L^x) \times [0, L^y] \times [-L^z, 0]``, (periodic, wall-bounded, wall-bounded) with values

|   Parameter             | Value       | Units | Description |
|   :-------:             | :---:       | :---:  |:---:       |
| ``L^x``           | ``10^3`` |  ``\text{km}``      | Zonal Length|
| ``L^y``           | ``10^{3}`` | ``\text{km}``           | Meriodional Length |
| ``L^z``           | ``3`` |  ``\text{km}``      | Depth |

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
\mathcal{S}^b(b,y,z) = - \frac{1}{\lambda^t}\left(b - \Delta b  \frac{\exp \left( z/h \right) - \exp \left( -L/h \right)}{\exp \left( 1 \right) - \exp \left( -L/h \right)} \right) \exp( \frac{y - L^y}{\lambda^N})
```

and the surface relaxation for buoyancy is

```math
\left. \Phi^*_b \right|_{\text{surface}} = \lambda^s \left(b - \Delta b \frac{y}{L^y} \right)
```

with parameter values

|   Parameter             | Value       | Units | Description |
|   :-------:             | :---:       | :---:  |:---:       |
| ``\lambda^s``           | ``10^{-4} `` |  ``\frac{m}{\text{s}}``      | surface buoyancy velocity constant|
| ``\Delta b``           | ``10 \times \alpha \times 10`` | ``\frac{\text{m}}{\text{ s}^2}``           | buoyancy jump|
| ``h``           | ``1`` | ``\text{km}``           | northern wall stratification e-folding length|
| ``\lambda^t``           | ``\frac{1}{28 \times 86400}`` | ``\frac{1}{\text{ s}}``           | relaxation time|
| ``\lambda^N``           | ``20`` | ``\text{km}``           | northern wall horizontal e-folding length|

The surface flux for the zonal velocity, (``u`` where ``\bm{u} = (u,v,w) ``) is
```math
\left. \Phi^*_u \right|_{\text{surface}} = -\frac{\tau}{\rho_0} \left( \exp \left( - \lambda^u \frac{(y-L/2)^2}{L^2} \right) - \exp \left( - \lambda^u \frac{(L/2)^2}{L^2} \right) \right)
```
with parameters

|   Parameter             | Value       | Units | Description |
|   :-------:             | :---:       | :---:  |:---:       |
| ``\tau``           | ``0.2 `` |  ``\frac{N}{\text{m}^2}``      | wind-stress magnitude|
| ``\lambda^u``           | ``32`` | dimensionless           | decay constant|

The purpose of the parameter ``\lambda^u`` is to remove any windstress being applied near the norther sponge relaxation. 

And finally, the linear drag was chosen to be
|   Parameter             | Value       | Units | Description |
|   :-------:             | :---:       | :---:  |:---:       |
| ``\mu``           | ``10^{-3}`` |  ``\frac{m}{\text{s}}``      | linear drag velocity|




