# Solar-Induced Fluorescence (SIF)

## Lee 2015 model

This page documents Lee 2015 model for the solar-induced fluorescence (SIF) module in ClimaLand, implemented in `src/standalone/Vegetation/solar_induced_fluorescence.jl`. SIF is computed at a single wavelength (755 nm) and represents the emission from the canopy surface.

## Summary
- **Inputs:** absorbed photosynthetically active radiation (APAR), photochemical yield (`Φ_PSII`), and parameters describing non‑photochemical quenching (NPQ).
- **Core:** computes fluorescence quantum yield from the balance between radiative decay, photochemistry, and non‑radiative processes.
- **Output:** canopy‑level SIF emitted at 755 nm (surface emission).

---

## Model formulation

### Fluorescence quantum yield
The fraction of excitations leading to fluorescence is

```math
Φ_F = \frac{k_F}{k_F + k_P + k_T},
```

where:
- `k_F` — radiative decay rate (fluorescence),
- `k_P` — photochemical quenching rate (charge separation in PSII),
- `k_T` — sum of constitutive and non‑photochemical quenching terms.

The photochemical yield is

```math
Φ_{PSII} = \frac{k_P}{k_F + k_P + k_T},
```

so that

```math
Φ_F = \frac{k_F}{k_F + k_T} (1 - Φ_{PSII}).
```

### Canopy emission
SIF at 755 nm is then calculated as

```math
SIF = Φ_F \, APAR,
```

where `APAR` is the absorbed PAR by the canopy (per ground area).
The module directly returns the **surface‑emitted fluorescence** without explicit partitioning into leaf‑level or top‑of‑canopy components.

---

## Variables & units

| Quantity | Symbol | Units | Notes |
|---|:--:|:--:|---|
| Fluorescence quantum yield | `Φ_F` | – | 0–1 |
| Photochemical yield (PSII) | `Φ_PSII` | – | input from photosynthesis model |
| Absorbed PAR | `APAR` | W m⁻² | absorbed by canopy |
| Canopy SIF | `SIF` | W m⁻² | emission at 755 nm |

---

## Parameters

| Parameter | Description | Typical value/range |
|---|---|---|
| `k_F` | Radiative decay rate | constant, relative scaling |
| `k_T` | Non‑radiative decay (constitutive + NPQ) | model parameter |
| `NPQ_params` | Non‑photochemical quenching parameters | model‑specific |

---

## Numerical details & coupling
- `Φ_PSII` is provided by the photosynthesis module (e.g., Farquhar‑type model).
- APAR is taken from canopy radiative transfer.
- NPQ can be specified via parameters or linked to dynamic formulations if available.
- The module returns the canopy surface SIF at 755 nm directly as a diagnostic.

---
