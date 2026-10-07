# Optimal LAI Model

The Optimal LAI model predicts seasonal to decadal dynamics of leaf area index (LAI) based on optimality principles, balancing energy and water constraints.

This model is based on [Zhou2025](@citet), which presents a general model for the seasonal to decadal dynamics of leaf area that combines predictions from both the light use efficiency (LUE) framework and optimization theory.

## Key Concepts and Definitions

### Potential vs Actual GPP

The model distinguishes between two types of gross primary productivity (GPP):

- **Potential GPP ($A_0$)**: The hypothetical GPP that would be achieved if the canopy absorbed all incoming photosynthetically active radiation (PAR). This corresponds to setting the fraction of absorbed PAR (fAPAR) to 1, which would require an infinitely dense canopy. This is computed as:
  ```math
  A_0 = \text{LUE} \times \text{PPFD}
  ```
  where LUE is the light use efficiency (mol CO₂ per mol photons) and PPFD is the photosynthetic photon flux density (mol photons m⁻² time⁻¹).

- **Actual GPP ($A$)**: The GPP achieved by the actual canopy with finite LAI:
  ```math
  A = A_0 \times \text{fAPAR} = A_0 \times (1 - e^{-k \cdot \text{LAI}})
  ```

### fAPAR and Beer-Lambert Law

The fraction of absorbed photosynthetically active radiation (fAPAR) follows Beer-Lambert's law:
```math
\text{fAPAR} = 1 - e^{-k \cdot \text{LAI}}
```
where $k$ is the light extinction coefficient. This represents the fraction of incoming PAR that is absorbed by the canopy.

## Model Overview

The optimal LAI model computes LAI dynamically by:
1. Calculating the seasonal maximum LAI (LAI$_{max}$) based on energy and water limitations
2. Computing steady-state LAI from daily meteorological conditions
3. Updating actual LAI using an exponential moving average to represent the lag in leaf development

## Seasonal Maximum LAI

The seasonal maximum LAI (LAI$_{max}$) is determined by the minimum of energy-limited and water-limited constraints (Equations 11-12 in Zhou et al. 2025):

```math
\begin{align}
\text{fAPAR}_{max} &= \min\left(\text{fAPAR}_{energy}, \text{fAPAR}_{water}\right) \\
\text{fAPAR}_{energy} &= 1 - \frac{z}{k \cdot A_{0,annual}} \\
\text{fAPAR}_{water} &= \frac{c_a(1-\chi)}{1.6 \cdot D_{growing}} \cdot \frac{f_0 \cdot P_{annual}}{A_{0,annual}} \\
\text{LAI}_{max} &= -\frac{1}{k} \ln(1 - \text{fAPAR}_{max})
\end{align}
```

**Physical interpretation:**
- **Energy-limited fAPAR**: Represents the optimal trade-off between carbon gain from photosynthesis and the cost of building/maintaining leaves. When $z / (k \cdot A_{0,annual})$ is large (high leaf cost relative to potential carbon gain), the optimal fAPAR is reduced.
- **Water-limited fAPAR**: Represents the constraint imposed by water availability. The numerator represents the water use efficiency (related to stomatal conductance), while the denominator relates to evaporative demand.

where:
- $A_{0,annual}$ is the annual total potential GPP (mol CO₂ m⁻² yr⁻¹) — the integrated daily $A_0$ over the year
- $P_{annual}$ is the annual total precipitation (mol H₂O m⁻² yr⁻¹). Conversion: 1 mm precipitation ≈ 55.5 mol H₂O m⁻²
- $D_{growing}$ is the mean vapor pressure deficit during the moist growing season (Pa): while the air is above 0 °C and the precipitation of the last 30 days is at least half its potential evaporation. Zhou et al. (2025) define the growing season by temperature alone; where it has a long dry season (seasonally dry tropics), the vegetation transpires mostly in the wet season, and the VPD of the dry season made the water limit too tight. Where the moist season is shorter than a month (deserts), $D_{growing}$ blends into the VPD of the whole season above freezing
- $k$ is the light extinction coefficient (dimensionless)
- $z$ is the unit cost of constructing and maintaining leaves (mol CO₂ m⁻² yr⁻¹): the costs of tree and grass leaves averaged geometrically with the tree share $t$ of the vegetation, $z = z_{tree}^{t} z_{grass}^{1-t}$. Grasses have the higher cost, as $z$ includes the below-ground allocation that supplies the leaves. The tree share is prescribed, by default from the natural vegetation of the CLM surface data, or computed from the simulated climate (`PrognosticTreeShare()`, see below)
- $c_a$ is the ambient CO₂ partial pressure (Pa). Conversion: 400 ppm at 101325 Pa ≈ 40 Pa
- $\chi$ is the ratio of leaf-internal to ambient CO₂ partial pressure (dimensionless), from stomatal optimization
- $f_0$ is the fraction of annual precipitation available to plants (dimensionless). It varies with the aridity index $AI$ as $f_0 = f_{0,max} \exp(-0.604 \ln^2(AI/1.9))$, peaking at $f_{0,max}$ at the energy–water transition

## Daily Steady-State LAI

Given daily meteorological conditions, the steady-state LAI ($L_s$) represents the LAI that would be in equilibrium with GPP if conditions were held constant (Equations 13-15):

```math
\begin{align}
\mu &= m \cdot A_{0,daily} \\
L_s &= \min\left\{\mu + \frac{1}{k} W_0[-k\mu \exp(-k\mu)], \text{LAI}_{max}\right\}
\end{align}
```

where:
- $A_{0,daily}$ is the daily potential GPP (mol CO₂ m⁻² day⁻¹)
- $m$ is a parameter relating steady-state LAI to steady-state GPP (Equation 20):

```math
m = \frac{\sigma \cdot \text{GSL} \cdot \text{LAI}_{max}}{A_{0,annual} \cdot \text{fAPAR}_{max}}
```

where GSL is the growing season length (days) and $\sigma$ is a dimensionless parameter representing departure from square-wave LAI dynamics (σ = 1 would mean LAI instantly reaches LAI$_{max}$ at the start of the growing season).

- $W_0$ is the principal branch of the Lambert W function, which satisfies $W(x) e^{W(x)} = x$

## LAI Update

The actual LAI is updated using an exponential weighted moving average to represent the time lag for photosynthate allocation to leaves (Equation 16):

```math
\text{LAI}_{new} = \alpha \cdot L_s + (1-\alpha) \cdot \text{LAI}_{prev}
```

where $\alpha$ is a smoothing factor (dimensionless, 0-1). The effective memory timescale is $\tau \approx 1/\alpha$ days. Setting $\alpha = 0.067$ corresponds to approximately 15 days of memory.

This LAI follows the seasonal optimum, so it falls to zero when the potential GPP does (in winter). The leaf area index used by the canopy keeps, for the tree share $t$ of the vegetation, at least a fraction $r$ (`optimal_lai_tree_retention`) of LAI$_{max}$ through the unfavourable season, the evergreen and semi-evergreen part of tree canopies:

```math
\text{LAI}_{canopy} = \text{LAI} + t \max(r \, \text{LAI}_{max} - \text{LAI}, 0)
```

with LAI$_{max}$ evaluated with the χ of the growing-season temperature and VPD. The evergreen share of trees is not modelled: it depends on biogeography and history more than on climate (spruce in boreal Canada, larch in eastern Siberia), and a share predicted from climate did not improve LAI over a constant one.

## Model Assumptions

1. **Water limitation enters once**: The potential GPP $A_0$ carries no soil moisture stress, so water availability acts only through the $f_0 P / A_0$ term of LAI$_{max}$ rather than being counted twice.
2. **Beer-Lambert light extinction**: Light absorption follows an exponential decay through the canopy.
3. **Optimal stomatal behavior**: The model assumes plants optimize their stomatal conductance following the P-model framework, giving the $\chi$ parameter.
4. **Growing season inputs are diagnosed from the simulated climate**: the growing season length is a trailing-year count of days above freezing, the growing-season VPD is the mean VPD while the air is above freezing, and $f_0$ follows the aridity index $AI = PET/P$ with $PET$ the Priestley-Taylor potential evapotranspiration of SPLASH (Davis et al. 2017), as in Zhou et al. (2025). Each is carried by a time-integrated variable in the prognostic state, seeded from the final state of a spin-up of the model (`experiments/long_runs/optimal_lai_spinup.jl`) or, where it has no data, from a climatology.
5. **Continuous update**: LAI and the trailing climate totals are advanced every timestep by the time-stepper, not by a daily callback.

## Parameters

| Parameter | Symbol | Unit | Typical Value | Description |
| :--- | :---: | :---: | :---: | :--- |
| Light extinction coefficient | $k$ | - | 0.5 | Controls light attenuation through canopy |
| Tree leaf cost | $z_{tree}$ | mol CO₂ m⁻² yr⁻¹ | 8.94 | Unit cost of building and maintaining tree leaves |
| Grass leaf cost | $z_{grass}$ | mol CO₂ m⁻² yr⁻¹ | 127 | Unit cost of building and maintaining grass leaves |
| LAI dynamics parameter | $\sigma$ | - | 1.08 | Departure from square-wave dynamics |
| Tree leaf retention | $r$ | - | 0.5 | Fraction of LAI$_{max}$ trees keep through the unfavourable season |
| Smoothing factor | $\alpha$ | - | 0.067 | Controls LAI response time (~15 days) |
| Peak precipitation fraction | $f_{0,max}$ | - | 0.65 | Fraction of precipitation used by plants at the energy–water transition |

The C3/C4 competition that sets the C3 fraction adds the coefficients below, fitted by
Lavergne et al. (2022) and used by pyrealm. The proportional C4 GPP advantage is passed
through a logistic, then penalised by the C3 tree cover $tc(g) = a g^b + c$ estimated from
the C3 GPP $g$ scaled to a year-long growing season (the annual total times 365/GSL, so
boreal forests keep the tree cover of their productivity during their short season), so
C4 is suppressed where C3 trees would shade it. With this
biomass model, the C3 fraction used by photosynthesis comes from the competition rather
than from the photosynthesis model's static map, so it requires the P-model. The
model computes its potential GPP and χ, and the per-pathway potential GPP the competition
compares, with its own P-model unit cost ratios, the pyrealm defaults with which Zhou et
al. (2025) and the competition were fitted, so that recalibrating the P-model's β for GPP
does not shift LAI or the C3/C4 balance.

| Parameter | Symbol | Unit | Typical Value | Description |
| :--- | :---: | :---: | :---: | :--- |
| Competition logistic steepness | $k_{c34}$ | - | 6.63 | Sharpness of the C4 fraction's response to the GPP advantage |
| Competition logistic midpoint | $q_{c34}$ | - | 0.16 | GPP advantage at which C3 and C4 are equally expected |
| Tree-cover coefficient | $a$ | - | 15.60 | Scale of the tree-cover relation |
| Tree-cover exponent | $b$ | - | 1.41 | Exponent of the tree-cover relation |
| Tree-cover offset | $c$ | - | -7.72 | Offset, so tree cover vanishes below a threshold GPP |
| Tree-cover reference GPP | $g_{ref}$ | kg C m⁻² yr⁻¹ | 2.8 | Normalizes the tree-cover relation to a proportion in [0, 1] |
| Unit cost ratio, C3 / C4 | $\beta_{C3}$, $\beta_{C4}$ | - | 146, 16.2 | P-model β of the model's potential GPP |

## Climate Tree Share

With `tree_share = PrognosticTreeShare()`, the tree share $t$ behind the leaf cost $z$ is
computed from the simulated climate rather than prescribed by a map, so it can change with
the climate:

```math
t = \frac{1}{1 + \exp[-(b_0 + b_L L_{tree} + b_d n_{dry} + b_T T_{growing})]}
```

where $L_{tree}$ is the LAI$_{max}$ a C3 tree canopy would reach (at $z_{tree}$, with the
χ of the growing-season temperature and VPD), $n_{dry}$ is the number of dry months of the
growing season (days above 0 °C when the precipitation of the last 30 days is below half
its potential evaporation, in months), and $T_{growing}$ is the mean air temperature while
above freezing (°C). Trees need water through the dry season and enough productivity to pay for
their canopy; grasses take over where the dry season is long. The coefficients, with the
leaf costs, σ and the tree leaf retention, were calibrated in an offline emulator of the
model (226 natural-vegetation points of a global 8° grid, 2008 ERA5) against MODIS LAI and
the tree share of the natural vegetation in the CLM5 surface data. Nothing in $t$
depends on the simulated LAI, so there is no feedback through the leaf cost. In this mode,
$t$ is also the tree share of the C3/C4 competition. The days and degree-days above
freezing, and the days and VPD of the moist growing season, are summed by time-integrated
variables whose memory grows from a day at the start of a simulation to $\tau_{long}$, so
that they are averages of their whole history until then.

| Parameter | Symbol | Unit | Typical Value | Description |
| :--- | :---: | :---: | :---: | :--- |
| Intercept | $b_0$ | - | -0.564 | Logistic intercept |
| Tree LAI coefficient | $b_L$ | m⁻² m² | 0.0464 | Effect of the LAI$_{max}$ of a tree canopy |
| Dry-month coefficient | $b_d$ | month⁻¹ | -0.364 | Effect of the number of dry months of the growing season |
| Temperature coefficient | $b_T$ | °C⁻¹ | 0.0621 | Effect of the growing-season temperature |

## Drivers

| Driver | Symbol | Unit | Description |
| :--- | :---: | :---: | :--- |
| Daily potential GPP | $A_{0,daily}$ | mol CO₂ m⁻² day⁻¹ | GPP assuming fAPAR = 1, without soil-moisture stress |
| Annual potential GPP | $A_{0,annual}$ | mol CO₂ m⁻² yr⁻¹ | Yearly integral of $A_0$ |
| Annual precipitation | $P_{annual}$ | mol H₂O m⁻² yr⁻¹ | Total yearly precipitation (1 mm ≈ 55.5 mol m⁻²) |
| Growing season VPD | $D_{growing}$ | Pa | Mean VPD during the moist growing season (T > 0°C, 30-day P ≥ PET/2) |
| Growing season length | GSL | days | Length of continuous period with T > 0°C |
| CO₂ partial pressure | $c_a$ | Pa | Ambient CO₂ (400 ppm ≈ 40 Pa at sea level) |

## Output

| Output | Symbol | Unit | Typical Range |
| :--- | :---: | :---: | :---: |
| Leaf Area Index | LAI | m² m⁻² | 0-10 |

## Implementation Notes

### Integration with Biomass Model

The optimal LAI model is implemented as a `ZhouOptimalLAIModel`, a subtype of `AbstractBiomassModel`. LAI is stored in `p.canopy.biomass.area_index.leaf`, consistent with the `PrescribedBiomassModel` interface. Auxiliary variables (A0 accumulators, GSL, precip, VPD, f0) are stored under `p.canopy.biomass.*`.

### Potential GPP Calculation

The implementation computes potential GPP ($A_0$) directly from the P-model with fAPAR = 1, the model's own unit cost ratios β, and no soil-moisture stress:

```math
A_0 = \text{PPFD} \times \text{LUE}
```

This differs from the global analysis of Zhou et al. (2025), whose $A_0$ includes the soil-moisture penalty of Stocker et al. (2020). Here water availability enters only through the $f_0 P / A_0$ term of LAI$_{max}$: a stressed $A_0$ would lower $A_0$ in dry regions and so loosen that limit.

### Daily and Annual A₀ Accumulation

- **Daily A₀**: Accumulated every timestep by integrating instantaneous A₀. Finalized at local noon.
- **Annual A₀**: Accumulated from daily values. Reset every 365 days.

### Unit Conversions

- **Precipitation**: 1 mm water = 1 kg m⁻² = 55.5 mol H₂O m⁻² (using molar mass of water = 18 g/mol)
- **CO₂ partial pressure**: At standard pressure (101325 Pa), 400 ppm CO₂ ≈ 40.5 Pa
- **A₀ units**: The P-model computes LUE in kg C/mol photons, so A₀ is converted to mol CO₂ using M$_c$ = 0.0120107 kg/mol

## References

[Zhou2025](@cite)
