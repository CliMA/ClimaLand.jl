## Energy+Hydrology

In more complex situations, describing the flow of liquid water in soil
is not sufficient. For example, to understand frozen soils, one must also
model phase changes, and keep track of the soil temperature to determine
when the water in the soil will freeze or melt. The more complex
soil model tracks liquid water, frozen water, and the energy of the soil,
and is able to capture phase change as well, by augmenting Richards Equation
with an equation for ice and energy.

We have

```math
\frac{\partial \vartheta_l}{\partial t} = - \nabla \cdot [-K \nabla h] + F_T/ρ_l + S_l\\
\frac{d \theta_i}{d t} = - F_T/ρ_i + S_i\\
\frac{\partial \rho e_{\rm{int}}}{\partial t} = - \nabla \cdot [-\kappa \nabla T - \rho e_{\rm{int,l}} K \nabla h]+S_e

```
where:
- $ϑ_l$ is the augmented volumetric liquid fraction, $t$ is the time, $K$ is the hydraulic conductivity, computed from $ϑ_l$ given a retention curve and a permeability curve, $ψ$ is the pressure head, which is computed from $ϑ_l$ given a retention curve function, and $h = ψ + z$ is the hydraulic head,
- $θ_i$ is the volumetric ice fraction,
- $ρe_{\rm{int}}$ is the volumetric internal energy, $κ$ is the thermal conductivity, $T$ is the temperature, $ρe_{int,l}$ is the volumetric internal energy of the soil liquid water,
- $S_e$, $S_i$, $S_l$ represent other sources of energy, ice, and water,
- $F_T$ is a source term representing phase changes, with $ρ_l$ and $ρ_i$ the density of ice and liquid water.

In order to solve these equations, the functions $ψ(ϑ_l,θ_i)$, $K(θ_i, ϑ_l, T)$, $κ(θ_i, ϑ_l)$, and
$T(ρe_int,θ_i, ϑ_l)$  must be specified. This in turn requires defining the
saturated conductivity $K_{\rm{sat}}$, the porosity $ν$, the residual water
content $θ_r$, and the parameters mapping saturation ϑ_l to K and ψ, and
the volumetric specific heat of the soil.

Other sources include root extraction of water (and corresponding extraction
of energy), subsurface runoff, and sublimation of ice.

ClimaLand supports both the van Genuchten and Brooks and Corey
retenton curve/permeability curve pairs, which we refer to in places
as the hydrology closure model. For the thermal conductivity, we use the model
of Balland and Arp (2003).

Since the liquid water and energy  partial differential equations ire stiff,
an implicit timestepping scheme must be used to advance them in time.