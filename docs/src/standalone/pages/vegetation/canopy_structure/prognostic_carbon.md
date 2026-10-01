# Prognostic Carbon Model

The prognostic carbon model (`PrognosticCarbonModel`) carries the live vegetation carbon
in four prognostic pools and predicts biomass from the carbon the canopy fixes. It wraps
an LAI model (`PrescribedBiomassModel` or `ZhouOptimalLAIModel`), which still sets the
leaf, stem and root area indices: GPP and LAI are the same with or without the pools.

## Pools and fluxes

The pools, in kg C m⁻² of ground, are the non-structural carbon (sugar) that GPP enters,
and the leaf, stem and fine-root carbon it is allocated to:

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
seasonally around it. The leaf and stem allocation fractions are blended between C3 and C4
values by the C3 fraction of the canopy, and roots take ``f_{root} = 1 - f_{leaf} - f_{stem}``.

## Climate dependence

There are no plant functional types. Two climate means, carried as time-integrated
variables with a two-year memory, modify the stem pool:

- **Precipitation**: dry climates build little wood (Sankaran et al., 2005). The stem
  allocation fraction is multiplied by ``x^n/(1+x^n)`` with ``x = P_{annual}/P_{half}``,
  and the allocation it withholds goes to roots.
- **Temperature**: trees live longer in cold climates. Below a reference mean annual
  temperature, the stem turnover time is multiplied by ``q^{(T_{ref} - T_{annual})/10}``,
  up to tenfold.

## Soil carbon

In `LandModel` and `SoilCanopyModel`, the litter makes the soil organic carbon prognostic:
``dSOC/dt = I_{litter}(z) - S_m``, with ``S_m`` the microbial respiration. Leaf and stem
litter enter on a shallow exponential profile and root litter follows the root
distribution, each normalized so the soil receives exactly the litter the pools shed.
Vegetation and soil carbon together then change by ``\text{GPP} - R_a - S_m``. These models
include the coupling (`SoilCarbonLitterInput`) when the canopy carries the pools; without
them, SOC is held at its initial condition.

## Usage

```julia
lai_model = Canopy.ZhouOptimalLAIModel{FT}(domain, toml_dict)
biomass = Canopy.PrognosticCarbonModel{FT}(lai_model, toml_dict)
canopy = Canopy.CanopyModel{FT}(domain, forcing, toml_dict; biomass)
```

The constructor then selects `PoolBasedAutotrophicRespirationModel` for the canopy
respiration, and the `cveg` diagnostic reports the pools. The pools start empty: with stem
turnover times of decades, an equilibrium state requires a spin-up.

## Parameters

| Parameter | Symbol | Unit | Default |
| :--- | :---: | :---: | :---: |
| Construction efficiency | ``a`` | - | 0.7 |
| Leaf / stem allocation, C3 | ``f_{leaf}``, ``f_{stem}`` | - | 0.3, 0.4 |
| Leaf / stem allocation, C4 | ``f_{leaf}``, ``f_{stem}`` | - | 0.4, 0.05 |
| Leaf / fine-root turnover time | ``\tau_{leaf}``, ``\tau_{root}`` | yr | 1.5, 2 |
| Stem turnover time, C3 / C4 | ``\tau_{stem}`` | yr | 30, 1 |
| Sapwood / fine-root respiration rate | ``r_{stem}``, ``r_{root}`` | yr⁻¹ | 0.1, 0.5 |
| Sapwood saturation | ``C_{sap,half}`` | kg C m⁻² | 2 |
| Maintenance respiration ``Q_{10}``, ``T_{ref}`` | | -, K | 2, 298.15 |
| Target sugar fraction | ``c_{nsc}`` | - | 0.1 |
| Allocation timescale | ``\tau_{alloc}`` | days | 10 |
| Precipitation at half stem allocation | ``P_{half}`` | m yr⁻¹ | 0.8 |
| Stem turnover temperature factor, reference | ``q``, ``T_{ref}`` | -, K | 2, 283 |
| Litter e-folding depth | | m | 0.05 |

The TOML names start with `carbon_`; see `toml/default_parameters.toml`.
