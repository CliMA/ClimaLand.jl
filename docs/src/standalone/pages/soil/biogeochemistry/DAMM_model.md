# Soil Biogeochemistry: CO₂ Production and O₂ Consumption

This section describes the coupled soil biogeochemistry model implemented in ClimaLand, which simulates the production and diffusion of CO₂, consumption of O₂, and microbial respiration sourced by a fixed pool of soil organic carbon ($Y_{C_{som}}$). The model combines the Dual Arrhenius and Michaelis-Menten (DAMM) kinetics framework for microbial respiration with gas diffusion equations that account for soil structure, temperature, and moisture, and with Henry's-law partitioning of dissolved CO₂ and O₂ between the gas and aqueous phases of the pore space.

## Model Overview

The biogeochemistry model carries three prognostic fields that evolve in space (with soil depth) and time. The variables, listed in the units used internally, are:

- the total CO₂ stored per unit soil volume (gas phase plus dissolved), $Y_{CO_2}$ (kg C m⁻³ soil), with $Y_{CO_2} = \theta_{\rm eff}\, c_g$ where $c_g$ is the air-phase carbon mass concentration (kg C m⁻³ air) and $\theta_{\rm eff}$ is the effective porosity, defined as the sum of the gas-filled pore volume and the dissolved-equivalent volume of the liquid pore water (see Effective Porosity below);
- the total O₂ stored per unit soil volume (gas phase plus dissolved), $Y_{O_2}$ (kg O m⁻³ soil), with $Y_{O_2} = \theta_{\rm {eff,o2}}\, c_{g, O_2}$ where $c_{g, O_2}$ is the air-phase carbon mass concentration (kg $O_2$ m⁻³ air) and $\theta_{\rm {eff, o2}}$ is the effective porosity, as above;
- the soil organic carbon density, $Y_{C_{som}}$ (kg C m⁻³), held fixed at its initial value (no tendency, see below).

The model captures the fundamental processes of aerobic decomposition: microbes consume soluble carbon and oxygen to produce carbon dioxide. CO₂ and O₂ diffuse through the soil pore space, with diffusion rates controlled by soil structure (porosity, tortuosity), moisture, temperature, and pressure.

$Y_{C_{som}}$ is initialized from the SoilGrids organic carbon density product and is treated as a stationary spatial parameter: we assume that the input of organic matter (e.g., from litter and root turnover) is equal to the rate of decomposition, so the actual $Y_{C_{som}}$ is constant in time. The model thus captures geographic heterogeneity of carbon stocks but does not integrate carbon depletion on simulation timescales — a deliberate simplification appropriate for the multi-decade runs targeted by ClimaLand. In the future, we may model the rate of decomposition explicitly and track $Y_{C_{som}}$ over time; expressing $Y_{C_{som}}$ as a prognostic variable now allows that change to be made easily.

## Governing Equations

### CO₂ Transport and Production

The prognostic variable $Y_{CO_2}$ is the **total** carbon mass stored per unit soil volume, summing the gas phase (occupying the air-filled pore volume $\theta_a$ at concentration $c_g$) and the dissolved phase (held in the liquid pore water $\theta_l$ at concentration $c_{aq}$, related to $c_g$ by Henry's law as $c_{aq} = \beta\, c_g$). The total is therefore

```math
\begin{equation}
Y_{CO_2} = \theta_a\, c_g + \theta_l\, c_{aq} = (\theta_a + \beta\, \theta_l)\, c_g \equiv \theta_{\rm eff}\, c_g,
\end{equation}
```

where $\theta_{\rm eff} \equiv \theta_a + \beta\, \theta_l$ is the **effective porosity** — the equivalent volume that would hold $Y_{CO_2}$ at concentration $c_g$ if all the carbon were in the gas phase. The gas-phase concentration is recovered as $c_g = Y_{CO_2} / \theta_{\rm eff}$. Storing $Y_{CO_2}$ rather than $c_g$ as the prognostic keeps the formulation well-behaved as the soil saturates ($\theta_a \to 0$): dissolved storage keeps $\theta_{\rm eff} \geq \beta\, \theta_l > 0$ as long as there is any liquid water. Diffusion is driven by gas-phase gradients:

```math
\begin{equation}
\frac{\partial Y_{CO_2}}{\partial t} = \nabla \cdot \left[D \, \nabla c_g\right] + S_m,
\qquad c_g = Y_{CO_2} / \theta_{\rm eff}
\end{equation}
```

where:

- the effective diffusivity of CO₂ in soil, $D$ (m² s⁻¹);
- the microbial CO₂ production rate, $S_m$ (kg C m⁻³ s⁻¹).

Here $\theta_a = \nu - \theta_l - \theta_i$ is the volumetric air content, i.e. the fraction of the soil volume occupied by the air phase. Note that we do not currently model fluxes of dissolved CO₂ due to fluxes of the soil water; only gas-phase diffusion is represented.

### O₂ Transport and Consumption

Oxygen transport and consumption in soil is described by:

```math
\begin{equation}
\frac{\partial Y_{\theta_{O_2}}}{\partial t}
= \nabla \cdot \left[D_{O_2}\, \nabla c_{g, O_2}\right]
- \frac{M_{O_2}}{M_C}\, S_m
\end{equation}
```

where:

- the effective porosity for O₂, $\theta_{\rm eff,O_2} = \theta_a + \beta_{O_2}\, \theta_l$ (m³ m⁻³ soil) — same construction as the CO₂ effective porosity defined above, but with the O₂-specific Henry's-law factor. Because O₂ is much less soluble than CO₂ ($\beta_{O_2} \approx 0.032$ vs $\beta_{CO_2} \approx 0.84$ at 298 K), the dissolved-equivalent term contributes much less for O₂;
- the O₂ mass concentration in air, $c_{O_2}$ (kg O₂ m⁻³ air);
- the effective diffusivity of O₂ in soil, $D_{O_2}$ (m² s⁻¹);
- the molar mass of O₂, $M_{O_2}$ (0.032 kg mol⁻¹);
- the molar mass of carbon, $M_C$ (0.012 kg mol⁻¹).

The second term represents O₂ consumption by microbial respiration. The stoichiometry follows C + O₂ → CO₂: each 12 g of carbon respired consumes 32 g of O₂. The conversion factor $M_{O_2}/ M_C)$ translates the carbon mass source $S_m$ (kg C m⁻³ s⁻¹ soil) into an oxygen source.

### Soil Organic Carbon

The $Y_{C_{som}}$ pool is held fixed at its initial value under the assumption of steady state:

```math
\begin{equation}
\frac{\partial Y_{C_{som}}}{\partial t} = 0
\end{equation}
```

The pool $Y_{C_{som}}$ enters the DAMM substrate kinetics (below) but is not depleted by respiration. Its value is set from the SoilGrids organic carbon density at the start of the simulation, providing the spatial structure of carbon availability.

## Microbial Respiration: DAMM Model

The microbial source term $S_m$ is computed using the Dual Arrhenius and Michaelis-Menten (DAMM) kinetics model [Davidson2012](@citet):

```math
\begin{equation}
S_m = V_{\max} \cdot MM_{sx} \cdot MM_{O_2}
\end{equation}
```

where:

- the maximum potential rate of respiration (temperature-dependent), $V_{\max}$ (kg C m⁻³ s⁻¹);
- the substrate availability factor, $MM_{sx}$ (dimensionless, 0–1);
- the oxygen limitation factor, $MM_{O_2}$ (dimensionless, 0–1).

### Maximum Respiration Rate

The maximum respiration rate follows Arrhenius kinetics, written here in centered form:

```math
\begin{equation}
V_{\max} = V_{{\rm ref},sx} \, \exp\!\left[-\frac{E_{a,sx}}{R}\!\left(\frac{1}{T} - \frac{1}{T_{{\rm ref},sx}}\right)\right]
\end{equation}
```

where:

- the maximum respiration rate at the reference temperature, $V_{{\rm ref},sx}$ (kg C m⁻³ s⁻¹);
- the reference temperature for the centered Arrhenius form, $T_{{\rm ref},sx}$ (K);
- the activation energy, $E_{a,sx}$ (J mol⁻¹).

This is algebraically equivalent to the uncentered form $V_{\max} = \alpha_{sx}\, \exp[-E_{a,sx}/(RT)]$ with $\alpha_{sx} = V_{{\rm ref},sx}\, \exp[E_{a,sx}/(R\, T_{{\rm ref},sx})]$. The centered form is preferred for calibration because $V_{{\rm ref},sx}$ (the rate at $T_{{\rm ref},sx}$) and $E_{a,sx}$ (the temperature sensitivity) are nearly orthogonal under parameter estimation, while $\alpha_{sx}$ and $E_{a,sx}$ are strongly correlated.

### Substrate Availability

The concentration of soluble carbon substrate available to microbes depends on diffusion through soil water films:

```math
\begin{equation}
[S_x] = p_{sx} \cdot Y_{C_{som}} \cdot D_{liq} \cdot \theta_l^{3}
\end{equation}
```

where:

- the concentration of soluble substrate, $[S_x]$ (kg C m⁻³);
- the soluble fraction of $Y_{C_{som}}$, $p_{sx}$ (dimensionless);
- the soil organic carbon density, $Y_{C_{som}}$ (kg C m⁻³);
- the dimensionless diffusion coefficient of soluble carbon in water, $D_{liq}$;
- the volumetric liquid water content, $\theta_l$ (m³ m⁻³).

The cubic dependence on moisture reflects the strong constraint that water films place on substrate diffusion to microbial sites.

The substrate limitation factor follows Michaelis-Menten kinetics:

```math
\begin{equation}
MM_{sx} = \frac{[S_x]}{kM_{sx} + [S_x]}
\end{equation}
```

where $kM_{sx}$ is the substrate Michaelis constant (kg C m⁻³).

### Oxygen Availability

The oxygen availability for microbial respiration accounts for diffusion limitations in the partially saturated pore space using a Millington–Quirk tortuosity model:

```math
\begin{equation}
O_{2,{\rm avail}} = D_{oa} \cdot O_{2,f} \cdot \theta_a^{4/3}
\end{equation}
```

where:

- the dimensionless O₂ availability metric, $O_{2,{\rm avail}}$;
- the dimensionless O₂ diffusion coefficient used in the availability scaling, $D_{oa}$;
- the volumetric O₂ fraction in air, $O_{2,f}$ (dimensionless) is computed from the mass concentration via the ideal gas law;
- the Millington–Quirk tortuosity factor, $\theta_a^{4/3}$.

The exponent 4/3 captures how tortuous diffusion pathways limit gas transport in partially saturated porous media.

The oxygen limitation factor follows Michaelis-Menten kinetics:

```math
\begin{equation}
MM_{O_2} = \frac{O_{2,{\rm avail}}}{kM_{O_2} + O_{2,{\rm avail}}}
\end{equation}
```

where $kM_{O_2}$ is the dimensionless O₂ Michaelis constant.

## Air–Water Partitioning of Dissolved Gases

CO₂ and O₂ partition between the soil's gas phase and the aqueous phase of pore water. Treating only the gas phase becomes singular as the soil saturates ($\theta_a \to 0$); accounting for the dissolved phase keeps the equations well-posed and is physically more accurate at high water content, especially for CO₂.

### Henry's Law

The temperature-dependent Henry's-law solubility is computed using a van 't Hoff form [Sander2015](@citet):

```math
\begin{equation}
K_H(T) = K_H(T_{\rm ref}^{H})\, \exp\!\left[\frac{\partial \ln K_H}{\partial T}\!\left(\frac{1}{T} - \frac{1}{T_{\rm ref}^{H}}\right)\right]
\end{equation}
```

where:

- the Henry's-law constant at the Sander reference temperature, $K_H(T_{\rm ref}^{H})$ (mol m⁻³ Pa⁻¹);
- the temperature coefficient (a dimensional analogue of $-\Delta_{\rm sol} H/R$), $\partial \ln K_H / \partial T$ (K);
- the Sander reference temperature, $T_{\rm ref}^{H} = 298.15$ K.

The dimensionless Henry's-law factor $\beta$ is defined by:

```math
\begin{equation}
\beta(T) = K_H(T)\, R\, T
\end{equation}
```

and converts a gas-phase concentration to an equivalent dissolved concentration: $c_{aq} = \beta\, c_g$. At 298 K this gives $\beta_{CO_2} \approx 0.84$ (dissolved CO₂ storage is comparable to gas-phase storage) and $\beta_{O_2} \approx 0.032$ (dissolved O₂ storage is small but stabilizing).

### Effective Porosity

The effective porosity used to translate between $Y_{CO_2}$ and $c_g$ — and similarly for the O₂ tendency — is the sum of the gas-filled pore volume and the dissolved-equivalent volume of the liquid pore water:

```math
\begin{equation}
\theta_{\rm eff} = \max(\theta_a + \beta\, \theta_l,\ \theta_{\rm eff,min})
\end{equation}
```

where the gas-filled pore volume $\theta_a$ and the volumetric liquid water content $\theta_l$ are introduced below. The term $\beta\, \theta_l$ is the "air-equivalent" volume of dissolved gas held in pore water; $\theta_{\rm eff,min} = 10^{-4}$ is a numerical floor that prevents division by zero in extreme conditions where both $\theta_a$ and $\theta_l$ are vanishingly small (see Numerical Safeguards). The same construction with $\beta_{O_2}$ defines $\theta_{\rm eff,O_2}$ used in the O₂ tendency.

## Gas Diffusivity in Soil

### CO₂ and O₂ Diffusivity

Both CO₂ and O₂ diffusivity in soil are computed using the same parameterization [Ryan2018](@citet), with separate reference diffusivities for the two gases:

```math
\begin{equation}
D_{\rm gas} = D_{0,{\rm gas}} \cdot (2\theta_{a100}^{3} + 0.04\, \theta_{a100}) \cdot \left(\frac{\theta_a}{\theta_{a100}}\right)^{2 + 3/b}
\end{equation}
```

where:

- the effective diffusivity of the gas in soil, $D_{\rm gas}$ (m² s⁻¹), evaluated for either CO₂ ($D$) or O₂ ($D_{O_2}$);
- the volumetric air content, $\theta_a$ (m³ m⁻³);
- the air-filled porosity at a soil water potential of $-100$ cm H₂O, $\theta_{a100}$ (dimensionless);
- the pore-size distribution parameter, $b$ (dimensionless).

The free-air reference diffusivity is corrected for soil temperature and pressure:

```math
\begin{equation}
D_{0,{\rm gas}} = D_{{\rm ref},{\rm gas}} \left(\frac{T}{T_{\rm ref}}\right)^{1.75} \left(\frac{P_{\rm ref}}{P}\right)
\end{equation}
```

where:

- the gas-specific reference diffusivity at standard conditions, $D_{{\rm ref},{\rm gas}}$ (m² s⁻¹), evaluated as $D_{\rm ref}$ for CO₂ and $D_{{\rm ref},O_2}$ for O₂;
- the standard reference temperature, $T_{\rm ref} = 273$ K;
- the standard reference pressure, $P_{\rm ref} = 101325$ Pa.

The temperature exponent 1.75 reflects increased molecular kinetic energy, while the pressure dependence accounts for the change in gas number density with surface pressure.

### Volumetric Air Content

The volumetric air content uses the total water content (liquid plus ice), so frozen water displaces air space:

```math
\begin{equation}
\theta_a = \nu - \theta_w, \qquad \theta_w = \theta_l + \theta_i
\end{equation}
```

where:

- the total soil porosity, $\nu$ (m³ m⁻³);
- the volumetric ice content, $\theta_i$ (m³ m⁻³).

By contrast, the dissolved storage in $\theta_{\rm eff}$ uses only $\theta_l$ — gases dissolve in liquid water, not in ice.

## O₂ Concentration Conversions

### From Volumetric Fraction to Mass Concentration

The O₂ mass concentration in air is computed from the volumetric fraction using the ideal gas law:

```math
\begin{equation}
c_{O_2} = O_{2,f} \frac{P\, M_{O_2}}{R\, T}
\end{equation}
```

where the O₂ mass concentration in air, $c_{O_2}$ (kg O₂ m⁻³ air), is recovered from the dimensionless volumetric fraction $O_{2,f}$ together with the surface pressure $P$, the molar mass $M_{O_2} = 0.032$ kg mol⁻¹, the gas constant $R = 8.314$ J mol⁻¹ K⁻¹, and the soil temperature $T$.

### From Mass Concentration to Volumetric Fraction

The inverse transformation is:

```math
\begin{equation}
O_{2,f} = c_{O_2}\, \frac{R\, T}{P\, M_{O_2}}
\end{equation}
```

These conversions are needed because diffusion is most naturally expressed in terms of mass concentration gradients, while the prognostic variable $Y_{\theta_{O_2}}$ and the microbial kinetics use the volumetric fraction.

## Boundary Conditions

### CO₂ Boundary Conditions

At the top boundary, the diffusive flux drives the soil gas-phase concentration toward the atmospheric value:

```math
\begin{equation}
F_{CO_2,{\rm top}} = -D\, \frac{c_{g,{\rm atm}} - c_{g,{\rm top}}}{\Delta z_{\rm top}},
\qquad c_{g,{\rm atm}} = \chi_{CO_2}\, P\, M_C / (R\, T_{\rm top})
\end{equation}
```

where:

- the atmospheric CO₂ mole fraction (mol mol⁻¹), $\chi_{CO_2}$, supplied as a driver;
- the gas-phase carbon concentration in the top soil layer, $c_{g,{\rm top}}$ (kg C m⁻³ air);
- the distance from the top layer center to the surface, $\Delta z_{\rm top}$ (m);
- the soil temperature at the top of the column, $T_{\rm top}$ (K).

At the bottom boundary, a zero-flux condition is applied: $F_{CO_2,{\rm bottom}} = 0$.

### O₂ Boundary Conditions

At the top boundary, the atmospheric O₂ volumetric fraction $\theta_{O_2,{\rm atm}}$ ($\approx 0.21$) is converted to a mass concentration via the ideal gas law and used as the boundary value:

```math
\begin{equation}
F_{O_2,{\rm top}} = -D_{O_2}\, \frac{c_{O_2,{\rm atm}} - c_{O_2,{\rm top}}}{\Delta z_{\rm top}},
\qquad c_{O_2,{\rm atm}} = \theta_{O_2,{\rm atm}}\, P\, M_{O_2} / (R\, T_{\rm top}).
\end{equation}
```

At the bottom boundary, a zero-flux condition is applied: $F_{O_2,{\rm bottom}} = 0$.

## Summary Tables

### Model Outputs

| Output | Symbol | Unit | Description |
| :--- | :---: | :---: | :--- |
| Total CO₂ | $Y_{CO_2}$ | kg C m⁻³ soil | Gas plus dissolved, $\theta_{\rm eff} c_g$ |
| Air-phase CO₂ | $c_g$ | kg C m⁻³ air | Diagnosed as $Y_{CO_2}/\theta_{\rm eff}$ |
| O₂ volumetric fraction | $O_{2,f}$ | dimensionless | O₂ fraction in soil air |
| Soil organic carbon | $Y_{C_{som}}$ | kg C m⁻³ | Initialized from SoilGrids, fixed |
| Microbial respiration | $S_m$ | kg C m⁻³ s⁻¹ | CO₂ production rate |

### Drivers

| Driver | Symbol | Unit | Range | Description |
| :--- | :---: | :---: | :---: | :--- |
| Soil temperature | $T$ | K | 250–320 | Controls respiration and diffusion |
| Volumetric liquid water | $\theta_l$ | m³ m⁻³ | 0.0–0.6 | Substrate diffusion, dissolved storage |
| Volumetric ice | $\theta_i$ | m³ m⁻³ | 0.0–0.6 | Displaces air in $\theta_a$ |
| Atmospheric pressure | $P$ | Pa | 80000–105000 | Affects gas diffusion |
| Atmospheric CO₂ | $\chi_{CO_2}$ | mol mol⁻¹ | $\sim 4 \times 10^{-4}$ | Top boundary condition |

### Parameters

| Parameter | Symbol | Unit | Typical Value | Description |
| :--- | :---: | :---: | :---: | :--- |
| Soil porosity | $\nu$ | m³ m⁻³ | 0.3–0.6 | Total pore space |
| Reference respiration rate | $V_{{\rm ref},sx}$ | kg C m⁻³ s⁻¹ | $\sim 2 \times 10^{-7}$ | DAMM, at $T_{{\rm ref},sx}$ |
| DAMM reference temperature | $T_{{\rm ref},sx}$ | K | 288.15 | Centered Arrhenius reference |
| Activation energy | $E_{a,sx}$ | J mol⁻¹ | $\sim 4 \times 10^{4}$ | Temperature sensitivity |
| Substrate Michaelis constant | $kM_{sx}$ | kg C m⁻³ | 0.01–0.5 | Half-saturation for substrate |
| O₂ Michaelis constant | $kM_{O_2}$ | dimensionless | 0.001–0.01 | Half-saturation for O₂ |
| Soluble $Y_{C_{som}}$ fraction | $p_{sx}$ | dimensionless | 0.005–0.5 | Fraction of $Y_{C_{som}}$ available |
| Pore size distribution | $b$ | dimensionless | 2–12 | Soil texture parameter |
| Air content at $-100$ cm H₂O | $\theta_{a100}$ | dimensionless | 0.05–0.3 | Reference air content |
| CO₂ Henry constant at $T_{\rm ref}^{H}$ | $K_H^{CO_2}(T_{\rm ref}^{H})$ | mol m⁻³ Pa⁻¹ | $3.4 \times 10^{-4}$ | Sander (2015) |
| O₂ Henry constant at $T_{\rm ref}^{H}$ | $K_H^{O_2}(T_{\rm ref}^{H})$ | mol m⁻³ Pa⁻¹ | $1.3 \times 10^{-5}$ | Sander (2015) |
| CO₂ Henry temperature coefficient | $\partial \ln K_H^{CO_2}/\partial T$ | K | 2400 | van 't Hoff slope |
| O₂ Henry temperature coefficient | $\partial \ln K_H^{O_2}/\partial T$ | K | 1500 | van 't Hoff slope |

### Constants

| Constant | Symbol | Unit | Value | Description |
| :--- | :---: | :---: | :---: | :--- |
| CO₂ diffusion coefficient (STP) | $D_{\rm ref}$ | m² s⁻¹ | $1.39 \times 10^{-5}$ | Reference diffusivity, CO₂ in air |
| O₂ diffusion coefficient (STP) | $D_{{\rm ref},O_2}$ | m² s⁻¹ | $1.67 \times 10^{-5}$ | Reference diffusivity, O₂ in air |
| Soil C substrate diffusivity | $D_{liq}$ | dimensionless | 3.17 | Substrate diffusion |
| O₂ availability scaling | $D_{oa}$ | dimensionless | 1.67 | Used in O₂ availability metric |
| Universal gas constant | $R$ | J mol⁻¹ K⁻¹ | 8.314 | Ideal gas law |
| Molar mass of O₂ | $M_{O_2}$ | kg mol⁻¹ | 0.032 | Oxygen molecular weight |
| Molar mass of C | $M_C$ | kg mol⁻¹ | 0.012 | Carbon atomic weight |
| Reference temperature (Sander) | $T_{\rm ref}^{H}$ | K | 298.15 | Henry's-law reference |
| Reference temperature (Ryan) | $T_{\rm ref}$ | K | 273 | Diffusivity reference |
| Reference pressure | $P_{\rm ref}$ | Pa | 101325 | Standard conditions |
| Effective porosity floor | $\theta_{\rm eff,min}$ | m³ m⁻³ | $10^{-4}$ | Numerical floor |
