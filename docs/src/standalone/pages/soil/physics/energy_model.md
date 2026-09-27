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

Since the liquid water and energy partial differential equations are stiff,
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
temperature $T_{\rm{top}}$ and thermal conductivity $κ_{\rm{top}}$) by a
thermal resistance

```math
r = \frac{Δz_{\rm{top}}}{κ_{\rm{top}}},
```

where $Δz_{\rm{top}}$ is the distance between the surface and the top cell
center, and $r_{\rm{litter}}$ (`litter_thermal_resistance`) is a resistance
supplied by the canopy biomass model (`canopy.biomass.r_litter`) for the
litter, thatch, and standing dead material between the skin and the mineral
soil, zero for bare soil. The skin temperature satisfies the surface energy
balance

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

The resistance $r$ is what separates the surface from the top cell center in
the discretization: for the 5 cm top layer of the global grid and a dry soil
($κ_{\rm{top}} ≈ 0.2$ W/m/K), $r ≈ 0.1$ m² K/W, and the skin is several
kelvin warmer than the top cell at midday. For the 2 cm layers of site
simulations $r$ is a few hundredths of m² K/W, which is still not negligible:
at US-Var in the dry season, the skin is on average about 3 K warmer than the
top cell at the daily maximum, and the diurnal range of the soil temperature
at 2–8 cm depth is about 10% smaller than when the fluxes are evaluated at the
top cell temperature. For standalone
`EnergyHydrology` and `SoilSnowModel`, the skin absorbs the downwelling
radiation with the soil albedo and emissivity and $r_{\rm{litter}} = 0$; in
`SoilCanopyModel` and `LandModel`, the shortwave and longwave radiation
reaching the soil are those transmitted and emitted by the canopy, and
$r_{\rm{litter}}$ takes its parameter value (0.1 m² K/W by default). Such a
resistance damps the diurnal cycle of the top soil layer at the expense of a
larger diurnal cycle of the skin, because the litter that supplies it also has
heat capacity, which the skin formulation neglects; the litter layer below
removes this limitation. Conduction between the soil and a
snowpack or lake sediment uses the top cell temperature.
### Litter layer

With a `SlabLitter` surface layer (`EnergyHydrology(...; surface_layer =
SlabLitter{FT}(toml_dict, Δt))`), a slab of litter with thickness $d_l$, thermal
conductivity $κ_l$, and volumetric heat capacity $ρc_l$ lies between the skin
and the soil, with a prognostic temperature $T_l$ (`Y.soil.T_litter`). The
thickness follows the canopy above,

```math
d_l = \max(d_{\rm{PAI}} \, ⟨PAI⟩, d_{\min}), \qquad
\frac{d⟨PAI⟩}{dt} = \frac{LAI + SAI - ⟨PAI⟩}{τ_{\rm{PAI}}},
```

where $⟨PAI⟩$ (`Y.soil.PAI_mean`) is an exponentially weighted trailing mean of
the plant area index with memory $τ_{\rm{PAI}}$ (`litter_memory_timescale`, one
year by default): litter stock is litterfall times residence time, both of which
average over phenology, so the litter is thickest under productive canopies and
persists when leaves are shed, with the single scaling $d_{\rm{PAI}}$
(`litter_thickness_per_pai`) for all land cover. Without a canopy, $⟨PAI⟩$
decays to zero and the slab reduces to the skin scheme at the floor $d_{\min}$.
The trailing mean is seeded with the canopy's plant area index at the start of a
simulation. The energy balance of the slab is

```math
\begin{aligned}
SW_n + LW_n(T_{\rm{sfc}}) + H(T_{\rm{sfc}}) + L(T_{\rm{sfc}}) &= \frac{T_l - T_{\rm{sfc}}}{r_{\rm{top}}}, &
r_{\rm{top}} &= \frac{d_l}{2 κ_l},\\
ρc_l d_l \frac{dT_l}{dt} &= -F_{\rm{atm}} + \frac{T_{\rm{top}} - T_l}{r_{\rm{bot}}}, &
r_{\rm{bot}} &= \frac{d_l}{2 κ_l} + \frac{Δz_{\rm{top}}}{κ_{\rm{top}}},
\end{aligned}
```

where $F_{\rm{atm}}$ is the net upward atmospheric flux at the skin (the left
side of the first equation). The soil then receives the conduction
$(T_{\rm{top}} - T_l)/r_{\rm{bot}}$ plus the fluxes that enter it directly: the
energy of infiltrating water and, under snow, the conduction from the snowpack,
which exchanges with the top soil cell so that the litter under snow follows the
soil. The total energy of the column includes $ρc_l d_l (T_l - T_0)$. As
$d_l \to 0$, this reduces to the skin scheme for bare soil. The slab damps the
diurnal cycle of both the top soil layer and the radiating skin.

The litter relaxes on a time scale $ρc_l d_l / (∂F_{\rm{atm}}/∂T_l)$ of a few
hundred seconds, shorter than the time step, so $T_l$ is treated implicitly. The
skin balance is solved once per step at $T_l^n$ and the atmospheric flux is
linearized about it, $F_{\rm{atm}} ≈ F_{\rm{atm}}^n + Λ (T_l - T_l^n)$ with
$Λ = ∂F_{\rm{atm}}/∂T_l$ from the skin solve. The backward Euler litter equation
is then linear in $T_l$ given $T_{\rm{top}}$ and is solved in closed form in
every Newton iteration of the soil solve, so that the soil top flux
$(T_{\rm{top}} - T_l(T_{\rm{top}}))/r_{\rm{bot}}$ depends on $T_{\rm{top}}$
alone; its derivative enters the soil energy Jacobian at the top face, which
keeps the litter–soil coupling fully implicit and the energy budget of soil plus
litter closed to the convergence of the Newton solve. This is why the
`SlabLitter` needs to be passed the time step, and why `LandSimulation` requires
the `ARS111` tableau with at least two Newton iterations for it. The linearized
flux is the exchange of the column with the atmosphere over the step and enters
the energy bookkeeping, together with the energy $(T_l - T_0)\,dC_l/dt$ of the
litter mass gained or lost as $⟨PAI⟩$ changes. The turbulent flux diagnostics,
and the fluxes a coupled atmosphere receives, are those at $T_l^n$, so the
land–atmosphere exchange is conservative to within $Λ (T_l^{n+1} - T_l^n)$ per
step. Litter water storage and the moisture dependence of $κ_l$ and $ρc_l$ are
not represented.

When the top cell is frozen (it contains ice and is below the depressed
freezing temperature), the skin temperature is capped at the depressed
freezing temperature: energy that would warm the skin further melts ice
instead. The cap is released once the top cell reaches the freezing
temperature, so that trace ice in a warm cell does not pin the skin.
The soil still receives the atmospheric fluxes at the capped skin temperature;
their excess over the skin–top conduction heats the top cell, where it melts
ice.

### Surface humidity and evaporation

The surface specific humidity used for evaporation is a conductance-weighted
mean of the specific humidity $q_{\rm{src}}$ of the vapor source at the skin
and the specific humidity of the air,

```math
q_{\rm{sfc}} = \frac{g_{\rm{soil}} \, q_{\rm{src}} + g_h \, q_{\rm{air}}}{g_{\rm{soil}} + g_h},
```

where $g_h$ is the aerodynamic conductance for heat and $g_{\rm{soil}}$ is the
conductance of the dry soil layer that forms at the surface as it dries
([SwensonLawrence2014](@citet)),

```math
g_{\rm{soil}} = \frac{D_v τ_a}{d_{\rm{sl}}}, \quad
d_{\rm{sl}} = d_{\rm{ds}} \left(\frac{α S_c - S_l}{α S_c}\right)^p \; \text{for } S_l < α S_c, \quad
τ_a = \frac{θ_a^{5/2}}{ν}, \quad θ_a = ν - θ_r - θ_i.
```

Here $D_v$ is the diffusivity of water vapor in air, $τ_a$ the tortuosity factor
for diffusion through the air-filled pore space of the dry layer
([Shokri2008](@citet)), whose liquid water content is residual (so that $θ_a$
does not depend on the moisture of the soil below), $S_l$ the effective saturation of the liquid water at
the surface (extrapolated from the top two layers, within the ice-free pore
space $ν - θ_i$), $S_c$ the critical saturation of the retention curve, and
$d_{\rm{ds}}$, $α$, and $p$ parameters ($p = 1$ is the linear form of
[SwensonLawrence2014](@citet)). When the surface is wetter than $α S_c$, no dry
layer exists, $g_{\rm{soil}}$ is unbounded, and $q_{\rm{sfc}} = q_{\rm{src}}$.
Conversely, $g_{\rm{soil}}$ is multiplied by $S_l^2 / (S_l^2 + S_0^2)$ with
$S_0 = 0.01$, so that it vanishes smoothly as the mobile liquid water of the top
layer is exhausted ($S_l → 0$) and only the immobile residual water remains.

The vapor source is the saturation specific humidity at the skin lowered by the
Kelvin factor of the soil water, $q_{\rm{src}} = h_r \, q_{\rm{sat}}(T_{\rm{sfc}})$
with $h_r = e^{g ψ_{\rm{sfc}} M_w / (R T_{\rm{sfc}})}$, when the air is drier
than that. When $q_{\rm{air}}$ lies between $h_r q_{\rm{sat}}$ and
$q_{\rm{sat}}$, the air is subsaturated over free water, so no dew forms, and
the Kelvin effect alone does not draw vapor into the soil: $q_{\rm{src}} =
q_{\rm{air}}$ and the vapor flux vanishes (as in CLM5). Dew forms when
$q_{\rm{air}} > q_{\rm{sat}}$; it condenses at the surface, so the dry-layer
resistance does not apply. With a `SlabLitter` surface layer, the diffusive
resistance of the litter slab, $r_{\rm{vap},l} = c_{\rm{vap}} (d_l - d_{\min}) / D_v$
with $c_{\rm{vap}}$ (`litter_vapor_resistance_factor`), is added in series with
$1 / g_{\rm{soil}}$. When the skin is below the (depressed) freezing
temperature, sublimation is computed instead, with $q_{\rm{sfc}}$ weighted by
the ice fraction $β_{\rm{ice}} = (θ_i / ν)^4$.

Under a canopy, the turbulent exchange of the ground is reduced to

```math
g_{\rm{eff}} = W g_h + \frac{1 - W}{1/g_h + r'}, \quad
W = e^{-(\rm{LAI} + \rm{SAI})}, \quad
r' = \frac{1 + 0.5 \min(\max(Ri, 0), 10)}{C_s u_*},
```

following CLM5: a fraction $W$ of the ground (the canopy gap fraction)
exchanges directly with the atmosphere, and the rest through the under-canopy
resistance $r'$ between the ground and the canopy air, with the transfer
coefficient $C_s$ (`undercanopy_ground_transfer_coefficient`), the friction
velocity $u_*$ above the canopy, and a stability correction in the bulk
Richardson number $Ri$ of the canopy air space. Since SurfaceFluxes.jl
evaluates the fluxes with $g_h$, the soil passes it an interface temperature
$T_i = T_a + (g_{\rm{eff}} / g_h)(T_{\rm{sfc}} - T_a)$, with $T_a$ the air
temperature brought dry-adiabatically to the surface, and the analogous
humidity, such that the fluxes from the interface with $g_h$ equal those from
the skin with $g_{\rm{eff}}$. Bare soil has $W = 1$ and $g_{\rm{eff}} = g_h$.

### Time treatment

The surface fluxes and the skin temperature are evaluated once per time step
from the state at the beginning of the step and held fixed during the implicit
solve for $ϑ_l$ and $ρe_{\rm{int}}$; without a litter layer they do not
contribute to the Jacobian, and with one only the litter–soil conduction does.
The column test used to check this shows the diurnal cycle of the top soil layer
to be insensitive to the time step at the step sizes used in ClimaLand
simulations.
