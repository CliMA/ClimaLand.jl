# Solar-Induced Fluorescence (SIF)

## Lee 2015 model

This page documents Lee 2015 model for the solar-induced fluorescence (SIF) module in ClimaLand, implemented in `src/standalone/Vegetation/solar_induced_fluorescence.jl`. SIF is computed at a single wavelength (755 nm) and represents the emission from the canopy surface. See [Lee2015](@citet).
The ClimaLand code page for this model can be found [here](https://github.com/CliMA/ClimaLand.jl/blob/main/src/standalone/Vegetation/solar_induced_fluorescence.jl).

---

## Model formulation

SIF at 755 nm is computed following Lee et al. (2015).

First, we compute the dark-adapted heat-loss rate coefficient $k_d$ as

```math
k_d = \max\left(k_{d,p1}(T_c - T_{freeze}) + k_{d,p2},\; k_{d,min}\right).
```

Here, $k_{d,p1}$ and $k_{d,p2}$ are parameters for heat loss in dark-adapted conditions (Tol et al. 2014, unitless), $k_{d,min}$ is the minimum allowed value of $k_d$, $T_c$ is canopy temperature (K), and $T_{freeze}$ is the freezing temperature threshold (K).

Next, we define

```math
x = 1 - \frac{J}{J_{max}},
```

where $J$ is the electron transport rate (mol m⁻² s⁻¹) and $J_{max}$ is the maximum electron transport capacity (mol m⁻² s⁻¹).

The light-adapted heat-loss rate coefficient $k_n$ is then given by

```math
k_n = (k_{n,p1}x - k_{n,p2})x,
```

where $k_{n,p1}$ and $k_{n,p2}$ are parameters for light-adapted conditions (Lee et al. 2013, unitless).

The baseline photochemical yield $\phi_{p0}$ is computed as

```math
\phi_{p0} = \frac{k_p}{\max(k_f + k_p + k_n,\; \varepsilon(FT))},
```

where $k_p$ is the rate coefficient for photochemical quenching, $k_f$ is the rate coefficient for fluorescence, and $\varepsilon(FT)$ denotes the machine epsilon for the floating-point type $FT$.

The effective photochemical yield is then

```math
\phi_p = \frac{J}{J_{max}}\,\phi_{p0}.
```

The fluorescence yield $\phi_f$ is computed as

```math
\phi_f = \frac{k_f}{\max(k_f + k_d + k_n,\; \varepsilon(FT))}(1 - \phi_p).
```

To convert from leaf-level fluorescence to observed signal, we compute

```math
\kappa = \kappa_{p1} V_{cmax25}^{leaf} \times 10^6 + \kappa_{p2},
```

where $V_{cmax25}^{leaf}$ is the maximum carboxylation rate at 25 °C (leaf level, mol m⁻² s⁻¹, internally converted to μmol), and $\kappa_{p1}, \kappa_{p2}$ are slope and intercept parameters from Lee et al. (2015).

The emitted fluorescence flux is

```math
F = APAR_{canopy}^{moles}\, \phi_f,
```

where $APAR_{canopy}^{moles}$ is the absorbed photosynthetically active radiation by the canopy (mol m⁻² s⁻¹).

Finally, the solar-induced fluorescence at 755 nm is

```math
SIF_{755} = \frac{F}{\max(\kappa,\; \varepsilon(FT))},
```

where $F$ is the emitted flux and $\kappa$ is the conversion factor relating leaf-level to observed fluorescence. $SIF_{755}$ is expressed in W m⁻².


## Model Parameters

| Symbol | Description | Units | Value |
|--------|-------------|-------|-------|
| $kf$ | Rate coefficient for fluorescence | unitless | 0.05 |
| $kd_{p1}$ | Parameter for heat loss in dark-adapted conditions (Tol et al. 2014) | unitless | 0.03 |
| $kd_{p2}$ | Parameter for heat loss in dark-adapted conditions (Tol et al. 2014) | unitless | 0.0273 |
| $min_{kd}$ | Minimum heat-loss coefficient in dark-adapted conditions (Tol et al. 2014) | unitless | 0.087 |
| $kn_{p1}$ | Parameter for heat loss in light-adapted conditions (Lee et al. 2013) | unitless | 6.2473 |
| $kn_{p2}$ | Parameter for heat loss in light-adapted conditions (Lee et al. 2013) | unitless | 0.5944 |
| $kp$ | Rate coefficient for photochemical quenching | unitless | 4.0 |
| $\kappa_{p1}$ | Slope relating leaf-level to observed fluorescence (Lee et al. 2015) | μmol⁻¹ m² s | 0.045 |
| $\kappa_{p2}$ | Intercept relating leaf-level to observed fluorescence (Lee et al. 2015) | unitless | 7.85 |