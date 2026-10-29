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
\frac{dC_{sugar}}{dt} &= \text{GPP} - R_m - S \\
\frac{dC_{leaf}}{dt}  &= a f_{leaf} S - C_{leaf}/\tau_{leaf} \\
\frac{dC_{stem}}{dt}  &= a f_{stem} S - C_{stem}/\tau_{stem} \\
\frac{dC_{root}}{dt}  &= a f_{root} S - C_{root}/\tau_{root}
\end{align}
```

``S`` is the allocation to growth, of which a fraction ``a`` becomes structure and
``(1-a)S`` is growth respiration, and the turnover terms are the litter passed to the soil.
The autotrophic respiration of the canopy is that of the pools, ``R_a = R_m + (1-a)S``.

### Maintenance respiration

```math
R_m = g\left(\frac{C_{sugar}}{C_{sugar,ref}}\right)\left[R_d + Q_{10}^{(T - T_{ref})/10}\left(r_{stem} C_{sap} + r_{root} C_{root}\right)\right],
\qquad C_{sap} = \frac{C_{stem}}{1 + C_{stem}/C_{sap,half}}
```

The leaf term is the dark respiration ``R_d`` of the photosynthesis model, which has its
own temperature response. The sapwood and fine roots respire at rates scaled by a
``Q_{10}`` of the canopy temperature ``T``. Only sapwood respires, so ``C_{sap}`` saturates
as the stem grows into mostly dead heartwood. The ramp ``g(x) = x^n/(1+x^n)`` shuts
respiration down as the sugar pool empties, so the pool cannot become negative.

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
``T_{\tau}``, the wood turnover time is multiplied by ``q^{(T_{\tau} - T_{annual})/10}``,
up to tenfold, with ``T_{annual}`` the air temperature averaged over two years (a
time-integrated variable).

## Soil carbon

In `LandModel` and `SoilCanopyModel`, the litter makes the soil organic carbon prognostic:
``dSOC/dt = I_{litter}(z) - S_m``, with ``S_m`` the microbial respiration. Leaf and stem
litter enter on a shallow exponential profile and root litter follows the root
distribution, each normalized so the soil receives exactly the litter the pools shed.
Vegetation and soil carbon together then change by ``\text{GPP} - R_a - S_m``. These models
include the coupling (`SoilCarbonLitterInput`) when the canopy carries the pools; without
them, SOC is held at its initial condition.

## Equilibrium and spin-up

The pools follow GPP and do not feed back on it. Under constant forcing, or periodic
forcing averaged over its period, each structural pool is its allocation times its
turnover time, ``C_i = a f_i \tau_i S``, and the allocation is the GPP left after
maintenance respiration; as the sapwood saturates with the wood pool, this is a quadratic
in ``S`` (`Canopy.equilibrium_carbon_pools`).

The pools start empty, and wood takes decades to fill. `spin_up_carbon_pools!` replaces
that transient: it advances a simulation over a spin-up period, averaging the drivers of
the pools, then sets the leaf, wood and fine-root pools to their equilibrium under them.
The sugar pool, which follows the fluxes of the last days, settles over weeks after such
a jump, so a second spin-up over the following year leaves a state without transients.

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
| Wood turnover temperature factor, reference | ``q``, ``T_{\tau}`` | -, K | 1.7, 283 |
| Litter e-folding depth | | m | 0.05 |

The wood allocation, turnover time and its temperature factor, and the fine-root
respiration rate are calibrated against the live woody carbon of [Xu2021](@citet) and
the ecosystem respiration of the carbon leaderboard (CarbonTracker NEE plus GOSIF GPP):
the wood turnover time is an effective value, which includes the disturbances that keep
the observed biomass below its potential. The TOML names start with `carbon_`; see
`toml/default_parameters.toml`.
