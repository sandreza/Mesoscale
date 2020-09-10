# Mesoscale Channel Setup
The equations of motion are [here](@ref sec:eom). 

# Boundary Conditions / Source Terms

- Buoyancy : Linear relaxation to the surface, insulating boundary conditions elsewhere. A sponge layer relaxation to the northern wall.
- Momentum : Linear drag and no penetration on the bottom, gaussian wind-stress for the zonal velocity at the surface, stress-free boundary conditions for the meridional velocity and no penetration at the surface

# Table of model parameters and functional form

The default viscocity and diffusivity values for the simulation are

|   Parameter             | Value       | Units | Description |
|   :-------:             | :---:       | :---:  |:---:       |
| ``\kappa^h``           | ``10^{-5} `` |  ``\frac{\text{m}^2}{\text{s}}``      | horizontal diffusivity |
| ``\kappa^v``           | ``10^{-2}`` | ``\frac{\text{m}^2}{\text{s}}``           | vertical diffusivity |
| ``\nu^h``           | ``10^{-5} `` |  ``\frac{\text{m}^2}{\text{s}}``      | horizontal viscocity |
| ``\nu^v``           | ``10^{-2}`` |   ``\frac{\text{m}^2}{\text{s}}``         | vertical viscocity |

The physical parameter values are

|   Parameter             | Value       | Units | Description |
|   :-------:             | :---:       | :---:  |:---:       |
| ``\rho^0``           | ``1024`` |  ``\frac{\text{kg}}{\text{m}^3}``      | reference density|
| ``\alpha``           | ``10^{-3}`` | ``\frac{\text{m}}{\text{s}^2 \text{K}}``           | thermal expansion coefficient |
| ``g``           | ``9.8 `` |  ``\frac{\text{m}}{\text{s}^2}``      | gravitational acceleration |
| ``\nu^v``           | ``10^{-2}`` |   ``\frac{\text{m}^2}{\text{s}}``         | vertical viscocity |

The domain is ``[0, L^x) \times [0, L^y] \times [-L^z, 0]``, (periodic, wall-bounded, wall-bounded) with values
|   Parameter             | Value       | Units | Description |
|   :-------:             | :---:       | :---:  |:---:       |
| ``L^x``           | ``1024`` |  ``\frac{\text{kg}}{\text{m}^3}``      | reference density|
| ``L^y``           | ``10^{-3}`` | ``\frac{\text{m}}{\text{s}^2 \text{K}}``           | thermal expansion coefficient |
| ``L^z``           | ``9.8 `` |  ``\frac{\text{m}}{\text{s}^2}``      | gravitational acceleration |

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
\mathcal{S}^b(b,y,z) = -\left(b - \Delta b  \frac{\exp \left( z/h \right) - \exp \left( -L/h \right)}{\exp \left( 1 \right) - \exp \left( -L/h \right)} \right) \exp()
```

and the surface relaxation for buoyancy is
```math
c = b - y/L\^y
```
with parameter values
|   Parameter             | Value       | Units | Description |
|   :-------:             | :---:       | :---:  |:---:       |
| ``h``           | ``20,000 `` |  ``\frac{1}{\text{m}}``      | vertical e-folding scale|
| ``\Delta b``           | ``10 \times \alpha \times g`` | ``\frac{1}{\text{m}} \frac{1}{\text{ s}}``           | linear term|


The 

