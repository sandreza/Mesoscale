# [Equations of Motions](@id sec:eom)

The equations of motion are the boussinesq equations with a linear equation of state,

```math
\begin{aligned}
\partial_t \bm{u} + \nabla \cdot ( \bm{u} \otimes \bm{u} ) + \bm{f} \times \bm{u} &= \nu^h \Delta^h \bm{u} + \nu^v \partial_z \bm{u} - \frac{1}{\rho_0}\nabla p+ b \hat{z}  \\ 
\nabla \cdot \bm{u} &= 0 \\
\partial_t b + \nabla \cdot (\bm{u} b) &= \kappa^h \Delta^h b + \kappa^v \partial_z b + \mathcal{S}^b
\end{aligned}
```
where we use an anistropic diffusivity and viscocity tensor, and ``\mathcal{S}^b`` is a buoyancy source term that will either be zero or a relaxation to a fixed temperature profile on the norther wall.
