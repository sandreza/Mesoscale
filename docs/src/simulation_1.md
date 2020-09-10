# Mesoscale Channel Setup
The equations of motion are [here](@ref sec:eom). 

# Boundary Conditions / Source Terms

- Buoyancy : Linear relaxation to the surface, insulating boundary conditions elsewhere. A sponge layer relaxation to the northern wall.
- Momentum : Linear drag and no penetration on the bottom, gaussian wind-stress for the zonal velocity at the surface, stress-free boundary conditions for the meridional velocity and no penetration at the surface

# Table of model parameters

The default viscocity and diffusivity values for the simulation are

|   Parameter             | Value       | Units | Description |
|   :-------:             | :---:       | :---:  |:---:       |
| ``\kappa^h``           | ``10^{-5} `` |  ``\frac{\text{m}^2}{\text{s}}``      | horizontal diffusivity |
| ``\kappa^v``           | ``10^{-2}`` | ``\frac{\text{m}^2}{\text{s}}``           | vertical diffusivity |
| ``\nu^h``           | ``10^{-5} `` |  ``\frac{\text{m}^2}{\text{s}}``      | horizontal viscocity |
| ``\nu^v``           | ``10^{-2}`` |   ``\frac{\text{m}^2}{\text{s}}``         | vertical viscocity |

We use a ``\beta``-plane coriolis 

with ``f`` and ``\beta`` values
|   Parameter             | Value       | Units | Description |
|   :-------:             | :---:       | :---:  |:---:       |
| ``f``           | ``-10^{-4} `` |  ``\frac{1}{\text{s}}``      | constant |
| ``\beta``           | ``10^{-11}`` | ``\frac{1}{\text{m}} \frac{1}{\text{ s}}``           | beta |

The parameters for the relaxation to the northern wall is 


The 

