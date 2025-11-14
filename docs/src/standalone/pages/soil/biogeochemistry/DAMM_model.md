# Soil Biogeochemistry: CO₂ Production and O₂ Consumption

This section describes the coupled soil biogeochemistry model implemented in ClimaLand, which simulates the production and diffusion of CO₂, consumption of O₂, and decomposition of soil organic carbon (SOC) through microbial respiration. The model combines the Dual Arrhenius and Michaelis-Menten (DAMM) kinetics framework for microbial respiration with gas diffusion equations that account for soil structure, temperature, and moisture dependencies.

## Model Overview

The biogeochemistry model tracks three prognostic variables that evolve in space (with soil depth) and time:
- CO₂ concentration in soil air ($C$, kg C m⁻³)
- Volumetric O₂ fraction in soil air ($O_{2a}$, dimensionless)
- Soil organic carbon ($C_{som}$, kg C m⁻³)

The model captures the fundamental processes of aerobic decomposition: microbes consume soil organic carbon and oxygen to produce carbon dioxide. These gases diffuse through the soil pore space, with diffusion rates controlled by soil structure (porosity, tortuosity), moisture content, and temperature.

## Governing Equations

### CO₂ Transport and Production

The evolution of CO₂ concentration in the soil follows a reaction-diffusion equation:

```math
\begin{equation}
\frac{\partial C}{\partial t} = -\nabla \cdot \left[-D \nabla C\right] + S_m
\end{equation}
```

where $C$ is the CO₂ concentration (kg C m⁻³), $D$ is the effective diffusivity of CO₂ in soil (m² s⁻¹), and $S_m$ is the microbial CO₂ production rate (kg C m⁻³ s⁻¹). The diffusive flux is $-D\nabla C$, directed from regions of high to low concentration.

### O₂ Transport and Consumption

Oxygen transport and consumption in soil is described by:

```math
\begin{equation}
\frac{\partial O_{2a}}{\partial t} = \frac{R T}{\theta_a P M_{O_2}} \nabla \cdot \left[D_{O_2} \theta_a \nabla \rho_{O_2}\right] - \frac{R T}{M_C \theta_a P} S_m
\end{equation}
```

where:
- $O_{2a}$ is the volumetric O₂ fraction in air (dimensionless, ~0.21 in atmosphere)
- $\theta_a$ is the volumetric air content (m³ air m⁻³ soil)
- $\rho_{O_2}$ is the O₂ mass concentration in air (kg O₂ m⁻³ air)
- $D_{O_2}$ is the effective diffusivity of O₂ in soil (m² s⁻¹)
- $R$ is the universal gas constant (8.314 J mol⁻¹ K⁻¹)
- $T$ is soil temperature (K)
- $P$ is atmospheric pressure (Pa)
- $M_{O_2}$ is the molar mass of O₂ (0.032 kg mol⁻¹)
- $M_C$ is the molar mass of carbon (0.012 kg mol⁻¹)

The second term represents O₂ consumption by microbial respiration, following the stoichiometry C + O₂ → CO₂, where each 12 g of carbon respired consumes 32 g of oxygen.

### Soil Organic Carbon Dynamics

The soil organic carbon pool is depleted by microbial decomposition:

```math
\begin{equation}
\frac{\partial C_{som}}{\partial t} = -S_m
\end{equation}
```

where $C_{som}$ is the soil organic carbon content (kg C m⁻³). This equation enforces carbon mass conservation: the rate of SOC consumption equals the rate of CO₂ production.

## Microbial Respiration: DAMM Model

The microbial source term $S_m$ represents heterotrophic respiration and is computed using the Dual Arrhenius and Michaelis-Menten (DAMM) kinetics model [Davidson2012](@citet):

```math
\begin{equation}
S_m = V_{max} \cdot MM_{sx} \cdot MM_{O_2}
\end{equation}
```

where $V_{max}$ is the maximum potential rate of respiration (temperature-dependent), $MM_{sx}$ represents substrate availability (0-1, dimensionless), and $MM_{O_2}$ is the oxygen limitation factor (0-1, dimensionless).

### Maximum Respiration Rate

The maximum potential respiration rate follows Arrhenius kinetics:

```math
\begin{equation}
V_{max} = \alpha_{sx} \exp\left(\frac{-Ea_{sx}}{RT}\right)
\end{equation}
```

where $\alpha_{sx}$ is the pre-exponential factor (kg C m⁻³ s⁻¹), $Ea_{sx}$ is the activation energy (J mol⁻¹), $R$ is the gas constant, and $T$ is soil temperature (K). This formulation captures the exponential increase in microbial activity with temperature.

### Substrate Availability

The concentration of soluble carbon substrates available to microbes depends on diffusion through soil water films:

```math
\begin{equation}
[S_x] = p_{sx} \cdot C_{som} \cdot D_{liq} \cdot \theta_l^3
\end{equation}
```

where:
- $[S_x]$ is the concentration of soluble substrate (kg C m⁻³)
- $p_{sx}$ is the fraction of SOC that is soluble (dimensionless)
- $C_{som}$ is the total soil organic carbon (kg C m⁻³)
- $D_{liq}$ is the diffusion coefficient of soluble carbon (dimensionless)
- $\theta_l$ is the volumetric liquid water content (m³ m⁻³)

The cubic dependence on moisture ($\theta_l^3$) reflects the strong constraint that water films place on substrate diffusion to microbes.

The substrate limitation factor follows Michaelis-Menten kinetics:

```math
\begin{equation}
MM_{sx} = \frac{[S_x]}{kM_{sx} + [S_x]}
\end{equation}
```

where $kM_{sx}$ is the Michaelis constant for substrate (kg C m⁻³).

### Oxygen Availability

The oxygen availability for microbial respiration accounts for diffusion limitations in porous media using a tortuosity model:

```math
\begin{equation}
O_{2,avail} = D_{oa} \cdot O_{2a} \cdot \theta_a^{4/3}
\end{equation}
```

where:
- $O_{2,avail}$ is the dimensionless O₂ availability metric
- $D_{oa}$ is the oxygen diffusion coefficient in air (dimensionless)
- $O_{2a}$ is the volumetric O₂ fraction in air
- $\theta_a^{4/3}$ is the Millington-Quirk tortuosity factor

The exponent 4/3 captures how tortuous diffusion pathways limit gas transport in partially saturated porous media.

The oxygen limitation factor follows Michaelis-Menten kinetics:

```math
\begin{equation}
MM_{O_2} = \frac{O_{2,avail}}{kM_{O_2} + O_{2,avail}}
\end{equation}
```

where $kM_{O_2}$ is the Michaelis constant for oxygen (dimensionless).

## Gas Diffusivity in Soil

### CO₂ and O₂ Diffusivity

Both CO₂ and O₂ diffusivity in soil are computed using the same parameterization [Ryan2018](@citet), which accounts for temperature, pressure, and soil moisture effects:

```math
\begin{equation}
D = D_0 \cdot (2\theta_{a100}^3 + 0.04\theta_{a100}) \cdot \left(\frac{\theta_a}{\theta_{a100}}\right)^{2 + 3/b}
\end{equation}
```

where:
- $D$ is the effective diffusivity in soil (m² s⁻¹)
- $\theta_a$ is the volumetric air content (m³ m⁻³)
- $\theta_{a100}$ is the air-filled porosity at soil water potential of -100 cm H₂O (dimensionless)
- $b$ is the pore size distribution parameter (dimensionless)

The reference diffusivity $D_0$ depends on temperature and pressure:

```math
\begin{equation}
D_0 = D_{ref} \left(\frac{T}{T_{ref}}\right)^{1.75} \left(\frac{P_{ref}}{P}\right)
\end{equation}
```

where $D_{ref}$ is the diffusion coefficient at standard conditions (m² s⁻¹), $T_{ref}$ = 273 K, $P_{ref}$ = 101325 Pa, $T$ is soil temperature (K), and $P$ is atmospheric pressure (Pa). The temperature dependence (exponent 1.75) reflects increased molecular kinetic energy, while the pressure dependence accounts for changes in gas density.

### Volumetric Air Content

The volumetric air content is computed as:

```math
\begin{equation}
\theta_a = \max(\nu - \theta_l, 0)
\end{equation}
```

where $\nu$ is the total soil porosity (m³ m⁻³) and $\theta_l$ is the volumetric liquid water content (m³ m⁻³). This simple relationship reflects that air fills the pore space not occupied by water.

## O₂ Concentration Conversions

### From Volumetric Fraction to Mass Concentration

The O₂ mass concentration in air is computed from the volumetric fraction using the ideal gas law:

```math
\begin{equation}
\rho_{O_2} = O_{2a} \frac{P \cdot M_{O_2}}{R \cdot T}
\end{equation}
```

where $\rho_{O_2}$ is the O₂ mass concentration in air (kg O₂ m⁻³ air), $O_{2a}$ is the volumetric O₂ fraction, $P$ is pressure (Pa), $M_{O_2}$ is the molar mass of O₂ (0.032 kg mol⁻¹), $R$ is the gas constant (8.314 J mol⁻¹ K⁻¹), and $T$ is temperature (K).

### From Mass Concentration to Volumetric Fraction

The inverse transformation is:

```math
\begin{equation}
O_{2a} = \rho_{O_2} \frac{R \cdot T}{P \cdot M_{O_2}}
\end{equation}
```

These conversions are necessary because diffusion is most naturally expressed in terms of mass concentration gradients, while the prognostic variable ($O_{2a}$) and the microbial kinetics use volumetric fractions.

## Boundary Conditions

### CO₂ Boundary Conditions

At the top boundary (soil-atmosphere interface), the model uses a diffusive flux condition with atmospheric CO₂ concentration:

```math
\begin{equation}
F_{CO_2,top} = -D \frac{C_{atm} - C_{top}}{\Delta z_{top}}
\end{equation}
```

where $C_{atm}$ is the atmospheric CO₂ concentration (kg C m⁻³), $C_{top}$ is the CO₂ concentration at the top soil layer, and $\Delta z_{top}$ is the distance from the top layer center to the surface.

At the bottom boundary, a no-flux condition is typically applied: $F_{CO_2,bottom} = 0$.

### O₂ Boundary Conditions

At the top boundary, atmospheric O₂ (volumetric fraction of 0.21) is used:

```math
\begin{equation}
F_{O_2,top} = -D_{O_2} \theta_{a,top} \frac{\rho_{O_2,atm} - \rho_{O_2,top}}{\Delta z_{top}}
\end{equation}
```

where $\rho_{O_2,atm}$ is computed from $O_{2a} = 0.21$ using the ideal gas law.

At the bottom boundary, a no-flux condition is typically applied: $F_{O_2,bottom} = 0$.

## Summary Tables

### Model Outputs

| Output | Symbol | Unit | Description |
| :---         |     :---:      |    :---:      |     :---   |
| CO₂ concentration | $C$ | kg C m⁻³ | CO₂ in soil air |
| O₂ volumetric fraction | $O_{2a}$ | dimensionless | O₂ fraction in soil air |
| Soil organic carbon | $C_{som}$ | kg C m⁻³ | Total SOC pool |
| Microbial respiration | $S_m$ | kg C m⁻³ s⁻¹ | CO₂ production rate |

### Drivers

| Driver | Symbol | Unit | Range | Description |
| :---         |     :---:      |    :---:      |     :---:   |    :---   |
| Soil temperature | $T$ | K | 250-320 | Controls respiration and diffusion |
| Volumetric liquid water | $\theta_l$ | m³ m⁻³ | 0.0-0.6 | Controls substrate diffusion |
| Atmospheric pressure | $P$ | Pa | 80000-105000 | Affects gas diffusion |
| Atmospheric CO₂ | $C_{atm}$ | kg C m⁻³ | - | Top boundary condition |

### Parameters

| Parameter | Symbol | Unit | Typical Range | Description |
| :---         |     :---:      |    :---:      |     :---:   |    :---   |
| Soil porosity | $\nu$ | m³ m⁻³ | 0.3-0.6 | Total pore space |
| Pre-exponential factor | $\alpha_{sx}$ | kg C m⁻³ s⁻¹ | 100e3-300e3 | DAMM kinetics |
| Activation energy | $Ea_{sx}$ | J mol⁻¹ | 50e3-70e3 | Temperature sensitivity |
| Substrate Michaelis constant | $kM_{sx}$ | kg C m⁻³ | 1e-10-0.1 | Half-saturation for substrate |
| O₂ Michaelis constant | $kM_{O_2}$ | dimensionless | 1e-10-0.1 | Half-saturation for O₂ |
| Soluble SOC fraction | $p_{sx}$ | dimensionless | 0.005-0.5 | Fraction of SOC available |
| Pore size distribution | $b$ | dimensionless | 2-12 | Soil texture parameter |
| Air content at -100 cm H₂O | $\theta_{a100}$ | dimensionless | 0.05-0.3 | Reference air content |

### Constants

| Constant | Symbol | Unit | Value | Description |
| :---         |     :---:      |    :---:      |     :---:   |    :---   |
| CO₂ diffusion coefficient (STP) | $D_{ref}$ | m² s⁻¹ | 1.47e-5 | Reference diffusivity |
| Soil C substrate diffusivity | $D_{liq}$ | dimensionless | 3.17 | Substrate diffusion |
| O₂ diffusion coefficient | $D_{oa}$ | dimensionless | 1.67 | O₂ diffusion in air |
| Universal gas constant | $R$ | J mol⁻¹ K⁻¹ | 8.314 | Ideal gas law |
| Molar mass of O₂ | $M_{O_2}$ | kg mol⁻¹ | 0.032 | Oxygen molecular weight |
| Molar mass of C | $M_C$ | kg mol⁻¹ | 0.012 | Carbon atomic weight |
| Reference temperature | $T_{ref}$ | K | 273 | Standard conditions |
| Reference pressure | $P_{ref}$ | Pa | 101325 | Standard conditions |
