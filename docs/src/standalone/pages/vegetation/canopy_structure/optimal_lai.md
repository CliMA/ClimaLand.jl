# Optimal LAI Model

The Optimal LAI model predicts seasonal to decadal dynamics of leaf area index based on optimality principles, balancing energy and water constraints.

This model is based on [Zhou2025](@citet), which presents a general model for the seasonal to decadal dynamics of leaf area that combines predictions from both the light use efficiency (LUE) framework and optimization theory.

## Model Overview

The optimal LAI model computes LAI dynamically by:
1. Calculating the seasonal maximum LAI (LAI$_{max}$) based on energy and water limitations
2. Computing steady-state LAI from daily meteorological conditions
3. Updating actual LAI using an exponential moving average to represent the lag in leaf development

## Seasonal Maximum LAI

The seasonal maximum LAI (LAI$_{max}$) is determined by the minimum of energy-limited and water-limited constraints:

```math
\begin{align}
\text{fAPAR}_{max} &= \min\left(\text{fAPAR}_{energy}, \text{fAPAR}_{water}\right) \\
\text{fAPAR}_{energy} &= 1 - \frac{z}{k \cdot A_{0,annual}} \\
\text{fAPAR}_{water} &= \frac{c_a(1-\chi)}{1.6 \cdot D_{growing}} \cdot \frac{f_0 \cdot P_{annual}}{A_{0,annual}} \\
\text{LAI}_{max} &= -\frac{1}{k} \ln(1 - \text{fAPAR}_{max})
\end{align}
```

where:
- $A_{0,annual}$ is the annual total potential GPP (mol m⁻² yr⁻¹)
- $P_{annual}$ is the annual total precipitation (mol m⁻² yr⁻¹)
- $D_{growing}$ is the mean vapor pressure deficit during the growing season (Pa)
- $k$ is the light extinction coefficient (dimensionless)
- $z$ is the unit cost of constructing and maintaining leaves (mol m⁻² yr⁻¹)
- $c_a$ is the ambient CO₂ partial pressure (Pa)
- $\chi$ is the ratio of leaf-internal to ambient CO₂ partial pressure (dimensionless)
- $f_0$ is the fraction of annual precipitation used by plants (dimensionless)

## Daily Steady-State LAI

Given daily meteorological conditions, the steady-state LAI ($L_s$) represents the LAI that would be in equilibrium with GPP if conditions were held constant:

```math
\begin{align}
\mu &= m \cdot A_{0,daily} \\
L_s &= \min\left\{\mu + \frac{1}{k} W_0[-k\mu \exp(-k\mu)], \text{LAI}_{max}\right\}
\end{align}
```

where:
- $A_{0,daily}$ is the daily potential GPP (mol m⁻² day⁻¹)
- $m$ is a parameter relating steady-state LAI to steady-state GPP, computed as:

```math
m = \frac{\sigma \cdot \text{GSL} \cdot \text{LAI}_{max}}{A_{0,annual} \cdot \text{fAPAR}_{max}}
```

where GSL is the growing season length (days) and $\sigma$ is a dimensionless parameter representing departure from square-wave LAI dynamics.

- $W_0$ is the principal branch of the Lambert W function

## LAI Update

The actual LAI is updated using an exponential weighted moving average to represent the time lag for photosynthate allocation to leaves:

```math
\text{LAI}_{new} = \alpha \cdot L_s + (1-\alpha) \cdot \text{LAI}_{prev}
```

where $\alpha$ is a smoothing factor (dimensionless, 0-1). Setting $\alpha = 0.067$ corresponds to approximately 15 days of memory.

## Parameters

| Parameter | Symbol | Unit | Typical Value | Description |
| :--- | :---: | :---: | :---: | :--- |
| Light extinction coefficient | $k$ | - | 0.5 | Controls light attenuation through canopy |
| Leaf construction cost | $z$ | mol m⁻² yr⁻¹ | 12.227 | Unit cost of building and maintaining leaves |
| CO₂ concentration ratio | $\chi$ | - | 0.7 | Ratio of leaf-internal to ambient CO₂ |
| Precipitation fraction | $f_0$ | - | 0.62 | Fraction of precipitation used by plants |
| LAI dynamics parameter | $\sigma$ | - | 0.771 | Departure from square-wave dynamics |
| Smoothing factor | $\alpha$ | - | 0.067 | Controls LAI response time (~15 days) |

## Drivers

| Driver | Symbol | Unit | Description |
| :--- | :---: | :---: | :--- |
| Daily potential GPP | $A_{0,daily}$ | mol m⁻² day⁻¹ | GPP with fAPAR = 1 |
| Annual potential GPP | $A_{0,annual}$ | mol m⁻² yr⁻¹ | Yearly integral of $A_0$ |
| Annual precipitation | $P_{annual}$ | mol m⁻² yr⁻¹ | Total yearly precipitation |
| Growing season VPD | $D_{growing}$ | Pa | Mean VPD when T > 0°C |
| Growing season length | GSL | days | Length of continuous period T > 0°C |
| CO₂ partial pressure | $c_a$ | Pa | Ambient CO₂ concentration |

## Output

| Output | Symbol | Unit | Range |
| :--- | :---: | :---: | :---: |
| Leaf Area Index | LAI | m² m⁻² | 0-10 |

## References

[Zhou2025](@cite)
