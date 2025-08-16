## Richards Equation

Richards Equation is a single variable partial differential equation
describing how water is distributed and flows in soil. In three-dimensional
space, we have

```math
\frac{\partial \vartheta_l}{\partial t} = - \nabla \cdot [-K \nabla h]
```
where $ϑ_l$ is the augmented volumetric liquid fraction, $t$ is the time, $K$ is the hydraulic conductivity, computed from $ϑ_l$ given a retention curve and a permeability curve, $ψ$ is the pressure head, which is computed from $ϑ_l$ given a retention curve function, and $h = ψ + z$ is the hydraulic head.

In order to solve this equation, the functions $ψ(ϑ_l)$ and $K(ϑ_l)$ must be specified. This in turn requires defining the
saturated conductivity $K_{\rm{sat}}$, the porosity $ν$, the residual water
content $θ_r$, and the parameters mapping saturation ϑ_l to K and ψ.

ClimaLand supports both the van Genuchten and Brooks and Corey
retenton curve/permeability curve pairs, which we refer to in places
as the hydrology closure model.

Since this partial differential equation is stiff, an implicit timestepping
scheme must be used to advance it in time.