# Slab Lake Model

The slab lake model represents inland water bodies (lakes, reservoirs) as a fixed-depth, fully mixed water column that exchanges energy with the atmosphere at its surface and with the underlying sediment at its base. The lake evolves energy prognostically and can undergo freezing and thawing, while maintaining a constant depth. In ClimaLand, the slab lake replaces the soil surface at grid points identified as inland water by a user-provided mask.

## Model Overview

The slab lake is a single-layer energy balance model. At each timestep:
1. Surface fluxes (radiation, sensible heat, latent heat, precipitation enthalpy) drive the lake energy budget
2. Heat is exchanged with the underlying soil sediment through thermal conduction
3. Runoff adjusts to maintain constant lake depth (the fixed-depth constraint)
4. Lake temperature and liquid fraction are diagnosed from internal energy

The key simplification is that the lake is **fully mixed**: there is no vertical temperature stratification within the water column.

## Governing Equation

The prognostic variable is the lake internal energy per unit ground area, $U$ (J m$^{-2}$):

```math
\frac{dU}{dt} = -Q_{\rm sfc} - Q_{\rm runoff} + Q_{\rm sed}
```

where we have:
- the net surface energy flux into the atmosphere, $Q_{\rm sfc}$ (W m$^{-2}$)
- the energy carried away by runoff, $Q_{\rm runoff}$ (W m$^{-2}$)
- the heat flux from the lake to the underlying sediment, $Q_{\rm sed}$ (W m$^{-2}$)

## Surface Energy Flux

The surface energy balance includes radiative, turbulent, and precipitation terms. Since $U$ is defined per unit ground area, all fluxes are scaled by the lake area fraction $f_{\rm lake}$:

```math
Q_{\rm sfc} = f_{\rm lake} \left( R_n + \text{LHF} + \text{SHF} + P_{\rm liq} \, \rho_l \, c_{p,l} \, (T_{\rm air} - T_0) \right)
```

where $R_n$ is the net radiation (adjusted for lake albedo), LHF and SHF are the latent and sensible heat fluxes, and the last term is the enthalpy flux of liquid precipitation arriving at the air temperature $T_{\rm air}$.

## Lake-Sediment Heat Exchange

Heat is exchanged between the lake mixed layer and the top soil layer through a series-conductance model:

```math
Q_{\rm sed} = -G_{\rm eff} \left(T_{\rm lake} - T_{\rm soil}\right)
```

where the effective conductance $G_{\rm eff}$ combines the lake conductance $G$ with the soil thermal conductance:

```math
G_{\rm eff} = \frac{1}{\frac{1}{G} + \frac{\Delta z_{\rm soil}}{2 \, \kappa_{\rm soil}}}
```

- the lake conductance, $G$ (W m$^{-2}$ K$^{-1}$), a single tunable parameter
- the top soil layer thermal conductivity, $\kappa_{\rm soil}$ (W m$^{-1}$ K$^{-1}$)
- the top soil layer thickness, $\Delta z_{\rm soil}$ (m)
- the lake and top soil temperatures, $T_{\rm lake}$ and $T_{\rm soil}$ (K)

The conductance $G$ encapsulates the lake-side thermal conductance (combining lake depth, conductivity, and mixing effects) into a single parameter that can be calibrated.

## Fixed-Depth Constraint and Runoff

The lake maintains a constant depth. The net surface water flux is:

```math
F_{\rm water} = P_{\rm liq} + E_{\rm liq} + E_{\rm ice}
```

where $P_{\rm liq}$ is the liquid precipitation rate, $E_{\rm liq}$ is the liquid evaporation rate, and $E_{\rm ice}$ is the ice sublimation rate. To enforce constant depth, the runoff is defined as:

```math
R_{\rm lake} = -F_{\rm water}
```

Positive runoff corresponds to net precipitation (water drains away), while negative runoff corresponds to net evaporation (water must be supplied to maintain depth). The associated energy flux is:

```math
Q_{\rm runoff} = R_{\rm lake} \cdot \rho \widehat{e}(T_{\rm lake}, q_l)
```

where $\rho \widehat{e}$ is the volumetric internal energy (blending ice and liquid properties according to $q_l$) evaluated at the lake temperature. Both runoff (effluent) and supplied water are assumed to carry energy at the lake temperature $T_{\rm lake}$. For runoff this is physically motivated: precipitation equilibrates with the lake before draining away as effluent. For evaporation-driven supply (negative runoff), the same assumption is adopted for simplicity, though in reality the source water temperature may differ.

## Thermodynamics and Phase Change

The lake temperature and liquid fraction are diagnosed from the internal energy. The model distinguishes three regimes.

### Fully frozen ($q_l \leq 0$)

```math
T_{\rm lake} = T_0 + \frac{U + \rho_l \cdot d \cdot L_{f,0}}{\rho_l \cdot d \cdot c_{p,i}}
```

### Mixed phase ($0 < q_l < 1$)

```math
T_{\rm lake} = T_{\rm freeze}
```

### Fully liquid ($q_l \geq 1$)

```math
T_{\rm lake} = T_0 + \frac{U}{\rho_l \cdot d \cdot c_{p,l}}
```

where we have:
- the reference temperature, $T_0$ (K)
- the density of liquid water, $\rho_l$ (kg m$^{-3}$)
- the lake depth, $d$ (m)
- the specific latent heat of fusion, $L_{f,0}$ (J kg$^{-1}$)
- the specific heat capacities of ice and liquid water, $c_{p,i}$ and $c_{p,l}$ (J kg$^{-1}$ K$^{-1}$)
- the freezing temperature, $T_{\rm freeze}$ (K)

### Liquid Fraction

The liquid fraction $q_l$ is determined by linear interpolation between the fully frozen and fully liquid energy states:

```math
q_l = \frac{U - U_{\rm ice}}{U_{\rm liq} - U_{\rm ice}}
```

where the bounding energies are:

```math
\begin{align}
U_{\rm ice} &= \rho_l \cdot d \cdot c_{p,i} \cdot (T_{\rm freeze} - T_0) - \rho_l \cdot d \cdot L_{f,0} \\
U_{\rm liq} &= \rho_l \cdot d \cdot c_{p,l} \cdot (T_{\rm freeze} - T_0)
\end{align}
```

$U_{\rm ice}$ is the internal energy when the lake is fully frozen at $T_{\rm freeze}$, and $U_{\rm liq}$ is the energy when fully liquid at $T_{\rm freeze}$.

## Surface Albedo

The lake surface albedo varies with phase state, interpolating linearly between open-water and ice-covered values:

```math
\alpha_{\rm lake} = q_l \cdot \alpha_{\rm liquid} + (1 - q_l) \cdot \alpha_{\rm ice}
```

Both the PAR and NIR albedo bands use this same value. Default values are $\alpha_{\rm liquid} = 0.08$ and $\alpha_{\rm ice} = 0.6$.

## Boundary Conditions at the Sediment

At lake grid points, the soil beneath the lake receives only conductive heat from the lake-sediment interface. No water infiltrates into the sediment:

- **Heat flux (top boundary)**: $Q_{\rm sed}$ (the lake-sediment heat flux)
- **Water flux (top boundary)**: 0 (no infiltration)
- **Surface/subsurface runoff**: 0 (runoff is handled by the lake)

All direct soil-atmosphere fluxes are suppressed at lake points, since the lake surface energy balance handles all exchanges with the atmosphere.

## Snow-Lake Coupling

When integrated with the snow model in a multi-component model, snow can accumulate on the lake surface. Snow-on-lake coupling is handled through the same fraction-based blending as the soil–lake interface. Detailed snow–lake partitioning of surface fluxes (e.g., separate treatment of exposed vs. snow-covered lake fractions) is not yet implemented and is planned for future work.

## Model Assumptions

1. **Well-mixed column**: The lake has no vertical temperature gradient. This is a reasonable approximation for shallow or wind-mixed lakes.
2. **Fixed depth**: Lake volume does not change. Excess water is removed as runoff; water deficits are implicitly supplied.
3. **No horizontal transport**: Each lake grid point is independent.
4. **Phase-dependent albedo**: Albedo transitions smoothly between open water and ice.
5. **No vegetation**: Canopy processes are not active at lake points.
6. **Isolated sediment**: The soil beneath the lake has zero water and energy flux at the top boundary (except for the lake-sediment heat flux), so it is hydrologically isolated regardless of its saturation state.

## Parameters

| Parameter | Symbol | Unit | Default Value | Description |
| :--- | :---: | :---: | :---: | :--- |
| Lake depth | $d$ | m | 10 | Mixed-layer depth |
| Conductance | $G$ | W m$^{-2}$ K$^{-1}$ | 0.1 | Lake–sediment conductance |
| Liquid albedo | $\alpha_{\rm liquid}$ | - | 0.08 | Open-water surface albedo |
| Ice albedo | $\alpha_{\rm ice}$ | - | 0.6 | Frozen lake surface albedo |
| Emissivity | $\varepsilon$ | - | 0.97 | Longwave emissivity |
| Momentum roughness length | $z_{0m}$ | m | 0.001 | Aerodynamic roughness |
| Scalar roughness length | $z_{0b}$ | m | 0.0001 | Heat/moisture roughness |

## Prognostic Variables

| Variable | Symbol | Unit | Description |
| :--- | :---: | :---: | :--- |
| Lake internal energy | $U$ | J m$^{-2}$ | Energy per unit ground area |

## Diagnostic Variables

| Variable | Symbol | Unit | Description |
| :--- | :---: | :---: | :--- |
| Lake temperature | $T_{\rm lake}$ | K | Diagnosed from $U$ |
| Liquid fraction | $q_l$ | - | 0 = fully frozen, 1 = fully liquid |
| Surface energy flux | $Q_{\rm sfc}$ | W m$^{-2}$ | Net flux into atmosphere |
| Sediment heat flux | $Q_{\rm sed}$ | W m$^{-2}$ | Heat exchange with soil |
| Lake runoff | $R_{\rm lake}$ | m s$^{-1}$ | Water flux to maintain constant depth |

## Implementation Notes

### Standalone InlandWater Component

The slab lake is implemented as a standalone component in `src/standalone/InlandWater/InlandWater.jl`, following the same pattern as Snow and Canopy. It is added to `LandModel` via the `prognostic_land_components` infrastructure when an `inland_water_mask` is provided.

Key types:
- `SlabLakeParameters{FT}`: lake physical parameters, including a single `conductance` parameter (W m⁻² K⁻¹) for sediment heat flux
- `SlabLakeModel{FT}`: the model struct, containing parameters, domain, boundary conditions, and the inland water mask
- `AtmosDrivenLakeBC`: boundary condition for atmosphere-driven lake simulations

The sediment heat flux uses a series resistance model combining the lake conductance with the soil thermal conductivity of the top soil layer.

### Fraction-Based Blending

At the `LandModel` level, soil boundary conditions are blended with lake fluxes using the inland water fraction `f`:
- `soil_water_bc = (1-f) * soil_water_bc`
- `soil_heat_bc = (1-f) * soil_heat_bc + f * sediment_heat_flux`


### Inland Water Mask

The inland water mask is a user-provided NetCDF file containing a fractional inland water coverage field. Grid points with non-zero lake fraction have their soil surface fluxes blended with lake fluxes.
