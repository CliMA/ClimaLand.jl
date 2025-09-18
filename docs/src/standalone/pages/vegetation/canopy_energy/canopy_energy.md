# Canopy Energy Balance (Bigleaf)

This page documents the **canopy energy** module used in ClimaLand's vegetation standalone runs, consistent with `src/standalone/Vegetation/canopy_energy.jl`.

---

## Summary
- **Purpose:** solve the canopy energy balance and return the **canopy (bigleaf) temperature**.
- **Inputs:** meteorology (air temperature, pressure, humidity, wind), radiation (shortwave & longwave), canopy/aerodynamic conductances, and stomatal control on water vapor transfer.
- **Outputs:** canopy temperature

---

## Model formulation (bigleaf)

### Net radiation (canopy surface)
```math
R_n = SW_{abs} + LW_{abs} - \varepsilon\,\sigma\,T_c^{4},
```
with absorbed shortwave/longwave `SW_{abs}`, `LW_{abs}`, emissivity `\varepsilon`, Stefan–Boltzmann constant `\sigma`, and canopy (bigleaf) temperature `T_c` (K).

### Turbulent fluxes (bulk)
Sensible heat:
```math
H = \rho\,c_p\,g_H\,(T_c - T_a),
```
Latent heat (bulk water‑vapour pathway):
```math
LE = \lambda\,g_w\,\big(q^*(T_c) - q_a\big),
```
where `\rho` is air density, `c_p` heat capacity, `\lambda` latent heat of vaporization, `g_H` the conductance to heat, `g_w` the conductance to water vapour (including stomatal effects), `q^*(T_c)` saturation specific humidity at `T_c`, and `q_a` ambient specific humidity.

### Energy balance closure
```math
R_n = H + LE + S,
```
with optional storage/biochemical term `S` (often negligible at standard timesteps). The code solves for `T_c` (hence `H`, `LE`) using a nonlinear iteration or linearized update as implemented.

---

## Variables & units

| Quantity | Symbol | Units | Notes |
|---|:--:|:--:|---|
| Canopy temperature | `T_c` | K | bigleaf surface temperature |
| Air temperature | `T_a` | K | input |
| Net radiation | `R_n` | W m⁻² | per ground area |
| Sensible heat | `H` | W m⁻² | per ground area |
| Latent heat | `LE` | W m⁻² | per ground area |
| Conductance to heat | `g_H` | m s⁻¹ | aero + boundary layer |
| Conductance to water | `g_w` | m s⁻¹ | includes stomatal pathway |
| Emissivity | `\varepsilon` | – | ~0.98 typical |

---

## Parameters

- `\varepsilon` — longwave emissivity of canopy leaves
- `\sigma` — Stefan–Boltzmann constant (5.670374419×10⁻⁸ W m⁻² K⁻⁴)
- Options controlling **aerodynamic stability corrections** for `g_H`
- Floors/tolerances for the **solver** used to update `T_c`

---
