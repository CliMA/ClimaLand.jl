# Prognostic Carbon Model

The prognostic carbon model (`PrognosticCarbonModel`) carries the live vegetation carbon
in four prognostic pools and predicts biomass from the carbon the canopy fixes. It wraps
an LAI model (`PrescribedBiomassModel` or `ZhouOptimalLAIModel`), which still sets the
leaf, stem and root area indices: GPP and LAI are the same with or without the pools.

## Pools and fluxes

The pools, in kg C m⁻² of ground, are the non-structural carbon (sugar) that GPP enters,
and the leaf, wood (stem) and fine-root carbon it is allocated to:

```math
\begin{align}
\frac{dC_{sugar}}{dt} &= M_C\,\text{GPP} - R_m - S \\
\frac{dC_{leaf}}{dt}  &= a f_{leaf} S - C_{leaf}/\tau_{leaf} \\
\frac{dC_{stem}}{dt}  &= a f_{stem} S - C_{stem}/(s\,\tau_{stem}) \\
\frac{dC_{root}}{dt}  &= a f_{root} S - C_{root}/\tau_{root}
\end{align}
```

where we have:
- the sugar (non-structural), leaf, wood and fine-root carbon, ``C_{sugar}``, ``C_{leaf}``, ``C_{stem}`` and ``C_{root}`` (kg C m⁻²)
- the gross primary production of the canopy, GPP (mol CO₂ m⁻² s⁻¹), from the photosynthesis model
- the molar mass of carbon, ``M_C`` (kg C mol⁻¹), which converts GPP to kg C
- the maintenance respiration, ``R_m`` (kg C m⁻² s⁻¹)
- the allocation to growth, ``S`` (kg C m⁻² s⁻¹)
- the allocation fractions of leaves, wood and fine roots, ``f_{leaf}``, ``f_{stem}`` and ``f_{root}`` (-), which sum to one
- the factor lengthening the wood turnover time in cold climates, ``s`` (-)

The allocation becomes structure with efficiency ``a``, and the rest, ``(1-a)S``, is growth
respiration. The turnover terms are the litter passed to the soil,
``L_i = C_i/\tau_i`` (kg C m⁻² s⁻¹), with ``s\,\tau_{stem}`` for wood. The autotrophic
respiration of the canopy is that of the pools, ``R_a = R_m + (1-a)S``. The sections below
give ``R_m``, ``S``, the fractions ``f_i`` and ``s``; the parameters, such as ``a`` and the
turnover times ``\tau_i``, are in the table at the end of the page.

### Maintenance respiration

```math
R_m = g\left(\frac{C_{sugar}}{C_{sugar,ref}}\right)\left[M_C R_d + Q_{10}^{(T - T_{ref})/10}\left(r_{stem} C_{sap} + r_{root} C_{root}\right)\right],
\qquad C_{sap} = \frac{C_{stem}}{1 + C_{stem}/C_{sap,half}}
```

where we have:
- the dark respiration of the leaves, ``R_d`` (mol CO₂ m⁻² s⁻¹), from the photosynthesis model, which has its own temperature response
- the canopy temperature, ``T`` (K)
- the sapwood carbon, ``C_{sap}`` (kg C m⁻²), the living part of the wood
- the ramp ``g(x) = x^n/(1+x^n)`` (-)

The sapwood and the fine roots both respire at rates scaled by a ``Q_{10}`` of the canopy
temperature, so that the formulation is the same in a standalone canopy, which has no soil
temperature. Only sapwood respires, and it saturates as the stem grows
into mostly dead heartwood. The ramp shuts respiration down as the sugar pool empties, so
the pool cannot become negative.

### Allocation

Allocation draws the sugar pool toward a target proportional to the living biomass,

```math
S = \frac{C_{sugar}}{\tau_{alloc}}\, g\left(\frac{C_{sugar}}{c_{nsc}\,(C_{leaf} + C_{sap} + C_{root})}\right),
```

so growth slows when reserves fall below the target, and the sugar pool oscillates
seasonally around it.

There are no plant functional types: the allocation fractions are those of trees and of
grasses, weighted by the tree share ``f_{tree}`` of the vegetation, and only trees build
wood,

```math
f_{leaf} = f_{tree}\, f_{leaf,tree} + (1 - f_{tree})\, f_{leaf,grass}, \qquad
f_{stem} = f_{tree}\, f_{stem,tree}, \qquad
f_{root} = 1 - f_{leaf} - f_{stem}.
```

The tree share is that of the LAI model: with `ZhouOptimalLAIModel`, the one that sets
its leaf cost, by default the climate tree share, a function of the seasonal maximum LAI
a tree canopy would reach, of the number of dry months of the growing season and of its
temperature; with prescribed LAI, it is given, e.g. the tree share of the natural
vegetation of the CLM surface data (`clm_tree_share`). Where the dry season is long,
the vegetation is grass and builds little wood, which separates savannas from forests
under the same annual rainfall.

### Wood turnover

Trees live longer in cold climates. Below a reference mean annual temperature
``T_{\tau}``, the wood turnover time is multiplied by

```math
s = \min\left(q^{\max(T_{\tau} - T_{annual},\, 0)/10},\ 10\right),
```

and a factor ``q \le 1`` disables the scaling (``s = 1``). The mean annual temperature
``T_{annual}`` (K) is a running mean of the air temperature ``T_{air}`` (K) with a memory
of two years, a time-integrated variable:

```math
\frac{dT_{annual}}{dt} = \frac{T_{air} - T_{annual}}{\tau_{climate}}.
```

## Soil carbon

In `LandModel` and `SoilCanopyModel`, the litter makes the soil organic carbon prognostic:
``dSOC/dt = I_{litter}(z) - S_m``, with the microbial respiration ``S_m`` of the
[Soil Biogeochemistry](../../soil/biogeochemistry/DAMM_model.md) model. Leaf and stem
litter enter on a shallow exponential profile and root litter follows the root
distribution, each normalized by its integral over the soil column so the soil receives
exactly the litter the pools shed:

```math
I_{litter}(z) = \left(L_{leaf} + L_{stem}\right)\frac{p(z; d_{litter})}{\int p(z; d_{litter})\,dz}
+ L_{root}\,\frac{p(z; d_{root})}{\int p(z; d_{root})\,dz},
\qquad p(z; d) = \frac{e^{z/d}}{d},
```

where we have:
- the soil organic carbon, SOC (kg C m⁻³), and the microbial respiration, ``S_m`` (kg C m⁻³ s⁻¹)
- the litter input to the soil organic carbon, ``I_{litter}`` (kg C m⁻³ s⁻¹)
- the leaf, wood and fine-root litter of the pools, ``L_{leaf}``, ``L_{stem}`` and ``L_{root}`` (kg C m⁻² s⁻¹)
- the depth coordinate, ``z`` (m), negative below the surface, the integrals being over the soil column
- the litter e-folding depth, ``d_{litter}`` (m), and the rooting depth of the LAI model, ``d_{root}`` (m)

Vegetation and soil carbon together then change by
``M_C\,\text{GPP} - R_a - \int S_m\,dz``. These models include the coupling
(`SoilCarbonLitterInput`) when the canopy carries the pools; without them, SOC is held at
its initial condition.

## Equilibrium and spin-up

The pools follow GPP and do not feed back on it. Under constant forcing, or periodic
forcing averaged over its period, each structural pool is its allocation times its
turnover time,

```math
C_{leaf} = a f_{leaf} \tau_{leaf} S, \qquad C_{stem} = k S, \qquad C_{root} = a f_{root} \tau_{root} S,
\qquad k = a f_{stem}\, s\, \tau_{stem},
```

and the allocation is the GPP left after maintenance respiration,
``S = M_C(\text{GPP} - R_d) - f_T\,(r_{stem} C_{sap} + r_{root} C_{root})``, with ``f_T``
the mean of ``Q_{10}^{(T - T_{ref})/10}``. The limitation of respiration by an empty sugar
pool (``g \approx 1``) is neglected. As the sapwood saturates with the wood pool, ``S``
solves

```math
(1 + \rho)\, S + \frac{f_T\, r_{stem}\, k S}{1 + k S / C_{sap,half}} = M_C\,(\text{GPP} - R_d),
\qquad \rho = f_T\, r_{root}\, a f_{root} \tau_{root},
```

a quadratic in ``S`` (`Canopy.equilibrium_carbon_pools`).

The pools start empty, and wood takes decades to fill. `spin_up_carbon_pools!` replaces
that transient: it advances a simulation over a spin-up period, averaging the drivers of
the pools (GPP, ``R_d``, ``f_T`` at the canopy temperature, the air temperature and the
tree share), then sets the leaf, wood and fine-root pools to their equilibrium under them,
and ``T_{annual}`` to the mean air temperature, which is also the ``T_{annual}`` of the
equilibrium. The sugar pool, which follows the fluxes of the last days, keeps its state
and settles over weeks after such a jump, so a second spin-up over the following year
leaves a state without transients.

## Usage

```julia
lai_model = Canopy.ZhouOptimalLAIModel{FT}(surface_domain, toml_dict)
biomass = Canopy.PrognosticCarbonModel{FT}(lai_model, toml_dict)
land = LandModel{FT}(forcing, toml_dict, domain, Δt; biomass, prognostic_land_components)
simulation = LandSimulation(start_date, stop_date, Δt, land)
for years in 1:2
    ClimaLand.Simulations.spin_up_carbon_pools!(simulation, start_date + Year(years))
end
solve!(simulation)
```

The pools then take the tree share of the LAI model; around a prescribed LAI, pass one,
e.g. `tree_share = Canopy.clm_tree_share(surface_space)`. The canopy respiration is that
of the pools (`PoolBasedAutotrophicRespirationModel`), and the `cveg` diagnostic reports
their sum, which the vegetation leaderboard compares with the live woody carbon of
[Xu2021](@citet).

## Parameters

| Parameter | Symbol | Unit | Default |
| :--- | :---: | :---: | :---: |
| Construction efficiency | ``a`` | - | 0.7 |
| Leaf / wood allocation of trees | ``f_{leaf,tree}``, ``f_{stem,tree}`` | - | 0.3, 0.32 |
| Leaf allocation of grasses | ``f_{leaf,grass}`` | - | 0.4 |
| Leaf / fine-root turnover time | ``\tau_{leaf}``, ``\tau_{root}`` | yr | 1.5, 2 |
| Wood turnover time | ``\tau_{stem}`` | yr | 24 |
| Sapwood / fine-root respiration rate | ``r_{stem}``, ``r_{root}`` | yr⁻¹ | 0.1, 0.75 |
| Sapwood saturation | ``C_{sap,half}`` | kg C m⁻² | 2 |
| Maintenance respiration ``Q_{10}``, ``T_{ref}`` | | -, K | 2, 298.15 |
| Target sugar fraction | ``c_{nsc}`` | - | 0.1 |
| Allocation timescale | ``\tau_{alloc}`` | days | 10 |
| Ramp exponent | ``n`` | - | 3 |
| Sugar pool of the respiration shutdown | ``C_{sugar,ref}`` | kg C m⁻² | 0.001 |
| Wood turnover temperature factor, reference | ``q``, ``T_{\tau}`` | -, K | 1.7, 283 |
| Memory of the mean annual temperature | ``\tau_{climate}`` | yr | 2 |
| Litter e-folding depth | ``d_{litter}`` | m | 0.05 |

The wood allocation, turnover time and its temperature factor, and the fine-root
respiration rate are calibrated against the live woody carbon of [Xu2021](@citet) and
the ecosystem respiration of the carbon leaderboard (CarbonTracker NEE plus GOSIF GPP):
the wood turnover time is an effective value, which includes the disturbances that keep
the observed biomass below its potential. The TOML names start with `carbon_`; see
`toml/default_parameters.toml`.
