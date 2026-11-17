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

Since the liquid water and energy  partial differential equations are stiff,
an implicit timestepping scheme must be used to advance them in time.

## Surface boundary conditions

When the soil is driven by the atmosphere (`AtmosDrivenFluxBC`), the boundary
fluxes at the soil surface are computed from the atmospheric and radiative
forcing and the soil surface state. The water flux is the sum of infiltration
(precipitation minus surface runoff) and evaporation, and the energy flux
(positive upward) is

```math
F = R_n + H + L + F_{\rm{infil}},
```

with the net radiation

```math
R_n = -(1 - α) SW_d - ϵ (LW_d - σ T_{\rm{sfc}}^4),
```

the sensible and latent heat fluxes $H$ and $L$ from Monin–Obukhov similarity
theory (SurfaceFluxes.jl), and $F_{\rm{infil}}$ the internal energy carried by
infiltrating water. All of these depend on the soil surface temperature
$T_{\rm{sfc}}$ and the surface specific humidity $q_{\rm{sfc}}$.

### Skin temperature

The radiating and turbulent-exchange surface of the soil is treated as a skin
with zero heat capacity, connected to the center of the top soil layer (at
temperature $T_{\rm{top}}$ and thermal conductivity $κ_{\rm{top}}$) by the
conduction resistance

```math
r = \frac{Δz_{\rm{top}}}{κ_{\rm{top}}} + \frac{Δz_{\rm{lit}}}{κ_{\rm{lit}}},
```

where $Δz_{\rm{top}}$ is the distance between the surface and the top cell
center, and the second term is the resistance of the litter on the soil, of
thickness $Δz_{\rm{lit}}$ and thermal conductivity $κ_{\rm{lit}}$ (see
Litter below; zero for a bare soil). The skin temperature satisfies the
surface energy balance

```math
SW_n + LW_n(T_{\rm{sfc}}) + H(T_{\rm{sfc}}) + L(T_{\rm{sfc}}) = \frac{T_{\rm{top}} - T_{\rm{sfc}}}{r},
```

with all fluxes positive upward, and is found by Newton's method within the
Monin–Obukhov iterations, in the same way as the snow surface temperature.
The same solve yields the turbulent fluxes at $T_{\rm{sfc}}$, which are
stored with it (`p.soil.turbulent_fluxes`), so one Monin–Obukhov solve per
step gives both. The soil column receives the atmospheric fluxes $F$
evaluated at $T_{\rm{sfc}}$, so energy is conserved even when the balance
closes only to the tolerance of the Monin–Obukhov solve. The skin temperature
is solved for with the atmospheric state at the reference height in
`p.drivers`, whether prescribed or supplied by a coupler; with a coupled
atmosphere, the same solve also provides the momentum and buoyancy fluxes.

For a bare soil, the resistance $r$ is what separates the surface from the
top cell center in the discretization: for the 5 cm top layer of the global
grid and a dry soil
($κ_{\rm{top}} ≈ 0.2$ W/m/K), $r ≈ 0.1$ m² K/W, and the skin is several
kelvin warmer than the top cell at midday. For the 2 cm layers of site
simulations $r$ is a few hundredths of m² K/W, which is still not negligible:
at US-Var in the dry season, the skin is on average about 3 K warmer than the
top cell at the daily maximum, and the diurnal range of the soil temperature
at 2–8 cm depth is about 10% smaller than when the fluxes are evaluated at the
top cell temperature. For standalone
`EnergyHydrology` and `SoilSnowModel`, the skin absorbs the downwelling
radiation with the soil albedo and emissivity; in `SoilCanopyModel` and
`LandModel`, the shortwave and longwave radiation reaching the soil are those
transmitted and emitted by the canopy. Conduction between the soil and a
snowpack or lake sediment uses the top cell temperature and does not include
the litter.

When the top cell is frozen (it contains ice and is below the depressed
freezing temperature), the skin temperature is capped at the depressed
freezing temperature: energy that would warm the skin further melts ice
instead. The cap is released once the top cell reaches the freezing
temperature, so that trace ice in a warm cell does not pin the skin.
The soil still receives the atmospheric fluxes at the capped skin temperature;
their excess over the skin–top conduction heats the top cell, where it melts
ice.

### Litter

The dead leaves, needles and twigs on the forest floor, and the dead leaves
under grasses, conduct heat poorly and insulate the mineral soil from the
surface. In a pine forest, the litter reduces the daily amplitude of the
temperature at the mineral soil surface by 5 to 6.7 K, depending on the soil
moisture ([Andrade2010](@citet)). Land surface models that resolve the forest
floor treat the litter as a layer with its own temperature and water content
([OgeeBrunet2002, HaverdCuntz2010, Napoly2017](@cite)). At three forest
sites, the composite ISBA model, which has no litter, overestimated the
diurnal amplitudes of the soil temperature and the ground heat flux; the
canopy and litter
layers of ISBA-MEB reduced the root mean square error of the ground heat flux
from 47 to 10 W/m², mainly owing to the canopy shading and the low thermal
diffusivity of the litter ([Napoly2017](@citet)).

ClimaLand represents the litter by its thermal resistance only. Its heat
capacity is neglected (that of 3 cm of dry litter is about 3% of that of the
5 cm top soil layer of the global grid), so the skin is the top of the litter,
and the resistance $r$ above is that of the conduction between the litter and
the soil in ISBA-MEB ([Napoly2017](@citet), Eq. A10), with the skin
temperature in place of the litter temperature. The thermal conductivity of
the litter, $κ_{\rm{lit}} = 0.1$ W/m/K, is the dry value of
[Napoly2017](@citet) (their Eq. A13); its increase with the water held by the
litter is not represented. The litter adds no resistance to the transfer of
water vapor: the soil evaporation is computed as described below, at the skin
temperature (in CLM4.5, by contrast, the litter is a resistance on the vapor
path only; [SakaguchiZeng2009](@citet)).

The standalone soil model has no litter by default (`NoLitter()`); passing a
`SoilLitter` as `litter` to the `EnergyHydrology` constructor adds one. In
`LandModel` and `SoilCanopyModel`, the default soil has the litter of
`ClimaLand.vegetation_litter`: wherever the maximum MODIS leaf area index is
positive, the litter thickness is

```math
Δz_{\rm{lit}} = f_{\rm{tree}} \, Δz_{\rm{tree}} + (1 - f_{\rm{tree}}) \, Δz_{\rm{grass}},
```

where the tree share $f_{\rm{tree}}$ is the cover of the tree plant functional
types over that of all vegetated types in the CLM surface data, and there is
no litter elsewhere. Under trees, $Δz_{\rm{tree}} = 3$ cm is the default
thickness of [Napoly2017](@citet), within the 1–10 cm reported for their
forest sites. No measurement was found for the litter under grasses: with
$Δz_{\rm{grass}} = 1$ cm, its conductance $κ_{\rm{lit}} / Δz_{\rm{grass}} = 10$
W/m²/K equals the skin conductivity of all vegetated tiles in the ECMWF land
surface model ([ECMWF2023IFS](@cite), Table 8.2), an empirical conductance
between the skin and the top soil layer (3.3 W/m²/K under 3 cm of litter).
Where the vegetation is all trees and the 5 cm top layer is dry, the litter
raises $r$ from about 0.1 to about 0.4 m² K/W.

At 79 FLUXNET sites, the litter reduces the median ratio of the simulated to
the observed diurnal amplitude of the soil temperature at the 2–5 cm depths of
the sensors from 2.6 to 1.1 when the depths are counted from the top of the
mineral soil, and to 1.5 when they are counted from the top of the litter.

### Surface humidity and evaporation

The surface specific humidity used for evaporation is a conductance-weighted
mean of the (soil water potential adjusted) saturation specific humidity at
the skin and the specific humidity of the air,

```math
q_{\rm{sfc}} = \frac{g_{\rm{soil}} \, q_{\rm{sat}}(T_{\rm{sfc}}) \, e^{g ψ_{\rm{sfc}} M_w / (R T_{\rm{sfc}})} + g_h \, q_{\rm{air}}}{g_{\rm{soil}} + g_h},
```

where $g_h$ is the aerodynamic conductance for heat and $g_{\rm{soil}}$ is
the conductance of the dry soil layer that forms at the surface as it dries
([SwensonLawrence2014](@citet)),

```math
g_{\rm{soil}} = \frac{D_v τ_a}{d_{\rm{sl}}}, \quad
d_{\rm{sl}} = d_{\rm{ds}} \left(\frac{α S_c - S_l}{S_l}\right)^p \; \text{for } S_l < α S_c, \quad
τ_a = \frac{θ_a^{5/2}}{ν}, \quad θ_a = ν - θ_r - θ_i.
```

Here $D_v$ is the diffusivity of water vapor in air, $τ_a$ the tortuosity
factor for diffusion through the air-filled pore space
([Shokri2008](@citet)) of the dry layer, whose liquid water content is
residual (so that $θ_a$ does not depend on the moisture of the soil below, as
in [SwensonLawrence2014](@citet)), $S_l$ the effective saturation of the
liquid water at the surface (extrapolated from the top two layers), $S_c$ the
critical saturation of the retention curve, and $d_{\rm{ds}}$, $α$, and $p$
parameters. When the surface is wetter than $α S_c$, no dry layer exists,
$g_{\rm{soil}}$ is unbounded, and $q_{\rm{sfc}}$ is the saturation value. When
the skin is below the (depressed) freezing temperature, sublimation is
computed instead, with $q_{\rm{sfc}}$ weighted by the ice fraction
$β_{\rm{ice}} = (θ_i / ν)^4$.

### Time treatment

The surface fluxes and the skin temperature are evaluated once per time step
from the state at the beginning of the step, held fixed during the implicit
solve for $ϑ_l$ and $ρe_{\rm{int}}$, and do not contribute to the Jacobian.
The column test used to check this shows the diurnal cycle of the top
soil layer to be insensitive to the time step at the step sizes used in
ClimaLand simulations.
