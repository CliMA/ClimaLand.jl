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
- $D_{growing}$ is the mean vapor pressure deficit during the growing season (Pa), where growing season is defined as days with T > 0°C
- $k$ is the light extinction coefficient (dimensionless)
- $z$ is the unit cost of constructing and maintaining leaves (mol CO₂ m⁻² yr⁻¹)
- $c_a$ is the ambient CO₂ partial pressure (Pa). Conversion: 400 ppm at 101325 Pa ≈ 40 Pa
- $\chi$ is the ratio of leaf-internal to ambient CO₂ partial pressure (dimensionless), from stomatal optimization
- $f_0$ is the fraction of annual precipitation available to plants (dimensionless), varies with aridity

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

## Model Assumptions

1. **Water limitation through soil moisture stress**: Both daily and annual potential GPP ($A_0$) include soil moisture stress (β). This allows vegetation structure (LAI$_{max}$, via $A_{0,annual}$) to adapt to water availability on annual timescales, while daily LAI dynamics respond to shorter-term moisture variability. This enables response to climate change and flushing events (with ~1 year lag for structural adaptation).
2. **Beer-Lambert light extinction**: Light absorption follows an exponential decay through the canopy.
3. **Optimal stomatal behavior**: The model assumes plants optimize their stomatal conductance following the P-model framework, giving the $\chi$ parameter.
4. **Growing season inputs are provided, not diagnosed**: The model takes growing season length (GSL) and growing-season VPD as inputs; it does not currently diagnose season onset/offset from temperature. Defaults assume a 240-day GSL.
5. **Daily update at local noon**: LAI is updated once per day at local solar noon.

## Parameters

| Parameter | Symbol | Unit | Typical Value | Description |
| :--- | :---: | :---: | :---: | :--- |
| Light extinction coefficient | $k$ | - | 0.5 | Controls light attenuation through canopy |
| Leaf construction cost | $z$ | mol CO₂ m⁻² yr⁻¹ | 12.227 | Unit cost of building and maintaining leaves |
| CO₂ concentration ratio | $\chi$ | - | 0.7 | Ratio of leaf-internal to ambient CO₂ |
| Precipitation fraction | $f_0$ | - | 0.62 | Fraction of precipitation used by plants |
| LAI dynamics parameter | $\sigma$ | - | 0.771 | Departure from square-wave dynamics |
| Smoothing factor | $\alpha$ | - | 0.067 | Controls LAI response time (~15 days) |

## Drivers

| Driver | Symbol | Unit | Description |
| :--- | :---: | :---: | :--- |
| Daily potential GPP | $A_{0,daily}$ | mol CO₂ m⁻² day⁻¹ | GPP assuming fAPAR = 1 with actual β |
| Annual potential GPP | $A_{0,annual}$ | mol CO₂ m⁻² yr⁻¹ | Yearly integral of daily $A_0$ with actual β |
| Annual precipitation | $P_{annual}$ | mol H₂O m⁻² yr⁻¹ | Total yearly precipitation (1 mm ≈ 55.5 mol m⁻²) |
| Growing season VPD | $D_{growing}$ | Pa | Mean VPD during growing season (T > 0°C) |
| Growing season length | GSL | days | Length of continuous period with T > 0°C; default 240 days |
| CO₂ partial pressure | $c_a$ | Pa | Ambient CO₂ (400 ppm ≈ 40 Pa at sea level) |

## Output

| Output | Symbol | Unit | Typical Range |
| :--- | :---: | :---: | :---: |
| Leaf Area Index | LAI | m² m⁻² | 0-10 |

## Implementation Notes

### Potential GPP Calculation

The implementation computes potential GPP ($A_0$) directly from the P-model with fAPAR = 1 and actual β (soil moisture stress):

```math
A_0 = \text{PPFD} \times \text{LUE}
```

where:
- PPFD is the total photosynthetic photon flux density (mol photons m⁻² s⁻¹), computed from downwelling PAR
- LUE is the light use efficiency with actual β (soil moisture stress)

This differs from Zhou et al. (2025) who use β = 1 for potential GPP. Our implementation uses actual β to allow vegetation structure (LAI$_{max}$) to respond to water availability on annual timescales, enabling response to climate change and interannual variability (e.g., flushing events with ~1 year lag).

The P-model intermediate values (ϕ₀, Γ*, η*, K$_{mm}$, ξ, c$_i$, m$_j$, m') are computed using the same formulations as in `pmodel.jl`.

### Daily and Annual A₀ Accumulation

- **Daily A₀**: Accumulated every timestep by integrating instantaneous A₀. Finalized at local noon.
- **Annual A₀**: Accumulated from daily values. Reset on January 1 of each year.

This approach ensures proper integration of A₀ over time rather than relying on approximations based on instantaneous values.

### Auxiliary Variables

The model maintains the following state in `p.canopy.lai_model`:
- `LAI`: Current leaf area index (m² m⁻²)
- `A0_daily`: Daily potential GPP from previous day (mol CO₂ m⁻² day⁻¹), with actual β
- `A0_annual`: Annual potential GPP from previous year (mol CO₂ m⁻² yr⁻¹), with actual β
- `A0_daily_acc`: Accumulator for current day's potential GPP (with actual β)
- `A0_annual_acc`: Accumulator for current year's potential GPP (with actual β)

### Unit Conversions

- **Precipitation**: 1 mm water = 1 kg m⁻² = 55.5 mol H₂O m⁻² (using molar mass of water = 18 g/mol)
- **CO₂ partial pressure**: At standard pressure (101325 Pa), 400 ppm CO₂ ≈ 40.5 Pa
- **A₀ units**: The P-model computes LUE in kg C/mol photons, so A₀ is converted to mol CO₂ using M$_c$ = 0.0120107 kg/mol

## References

[Zhou2025](@cite)
