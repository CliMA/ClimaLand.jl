# Stomatal conductance
Stomata play an important role in uptaking CO2 for photosynthesis while limiting water loss during transpiration. Consequently, an accurate depiction of stomatal conductance is required to study leaf energy fluxes, transpiration, and photosynthesis.

This section describes multiple models of stomatal conductance implemented in ClimaLSM. 

## Medlyn Model

The Medlyn model is a semiempirical model that relates stomatal conductance and photosynthesis and is derived from water-use efficiency optimization theory.

Transpiration is computed using the stomatal conductance and Monin-Obukhov theory.
```math
\begin{equation}
T = -\rho_a g_{\mathrm{eff}} \left[q_{a}- q_v(T_\mathrm{leaf}, \rho_\mathrm{sfc}) \right],
\end{equation}
```
where $T$ is the transpiration (mass flux of water vapor), $q_{a}$ is the specific humidity at the lowest level of the atmosphere, $q_v(T_\mathrm{T_{leaf}}, \rho_\mathrm{sfc})$ is the saturated specific humidity over liquid water, given the temperature of the leave $T_{leaf}$ and air density at the surface $\rho_{sfc}$. We will approximate $T_{leaf} = T_{a}$ and $\rho_{\mathrm{sfc}} = \rho_a$.

We also need the effective conductivity, given by
```math
\begin{equation}
    g_{\mathrm{eff}} = \frac{1}{g_{\mathrm{ae}}^{-1}+g_{\mathrm{s}}^{-1}},
\end{equation}
```
where $g_{ae}$ is the aerodynamic conductance, computed by the MOST solve, and $g_s$ is the stomatal conductance to water vapor per unit ground area. The units of all conductances are $m/s$.

The stomatal conductance is calculated using the Medlyn stomatal conductance model (Medlyn, 2011), while omitting cuticular and epidermal losses by assuming zero minimum stomatal conductance:
```math
\begin{align}
g_{s,m}(PAR, T, VPD, c_a) &= g_{0,m} + D_{rel} \times m \frac{A_n(PAR, T, VPD, c_a)}{c_a}\nonumber \\
g_s &= \frac{g_{s,m}}{\rho_m}
\end{align}
```
where $D_{rel} =1.6$ (unitless) is the relative diffusivity of water vapor with respect to CO$_2$, $\rho_m$ is the molar density of water, and $m$ is the Medlyn factor,
```math
\begin{equation}
    m = \left( 1 + \frac{g_1}{\sqrt{VPD}} \right), 
\end{equation}
```
where g1 is the slope parameter, inversely proportional to the square root of marginal water use efficiency (Medlyn, 2011). We also have A$_n$ as the biochemical demand for CO$_2$ calculated using the photosynthesis model (Farquhar, 1980; Equation \eqref{eq:a_n}; units of molar flux). The resulting units are $m/s$. $g_{0,m}$ is a minimum molar conductivity. (subscript $m$ indicates molar).

The model has the following parameters:

| Constants | Symbol | Unit | Value |
| :---         |     :---:      |    :---:      |     :---:   |
| Relative diffusivity of water vapor | $D_{rel}$  | - | 1.6 |
| Minimum stomatal conductance | $g_0$ | mol/$m^2$/s | 1e-4 |
| Slope parameter | $g_1$  | $\sqrt{Pa}$  | 790 |

