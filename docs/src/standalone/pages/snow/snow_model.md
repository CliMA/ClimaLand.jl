## Snow modeling

The ClimaLand snow model is composed of a simple processed based model for conservation
of water mass and energy and parameterizations for albedo,
depth, and surface temperature.

The key conservative equations are

```math
\rho_l\frac{\partial S}{\partial t} =  P_{\rm {snow}} + \sigma(P_{\rm liq} - E - \rho_l R_{\rm snow}) \\
\frac{\partial U}{\partial t} & =  e_{i, \rm snow} P_{\rm {snow}} + \bigg( -\rho_l e_l R_{\mathrm{snow}} + e_{l, \rm rain} P_{\rm liq} +G -H - LE + R_n\bigg)\sigma \\
\rho_l \frac{\partial S_l}{\partial t} =  \bigg(P_{\rm liq} - E_l + \rho_lM - \rho_lR_{\rm snow}\bigg) \sigma,
```
where:
- $S$ is the snow water equivalent per ground area, $t$ is the time, $P_{\rm {snow}}$ and $\sigma(P_{\rm liq}$ are the precpitation fluxes, $E$ is the vapor flux, and $\rho_l R_{\rm snow}$ is a mass flux of water due to meltwater leaving the snowpack,
- $U$ is the energy of the snowpack per unit ground area, and $G$ is the snow/soil ground heat flux, $H$ is the sensible het flux, $LE$ is the latent heat flux, and $R_n$ is the radiative heat flux, with terms like $e_x$ reflecting the energy associated with precipitation or runoff, with $x$ indicating snowfall, rain, or runoff of liquid water.
- $S_l$ is the liquid water equivalent per ground area, and $\rho_lM$ a flux associated with phase changes within the snowpack,
- $\sigma$ is the snow cover fraction.

In order to solve these equations, we must define parameterizations for computing the fluxes ($H, LE, R_n, M, R$).
Further detail on these will be added as the model details are finalized. 