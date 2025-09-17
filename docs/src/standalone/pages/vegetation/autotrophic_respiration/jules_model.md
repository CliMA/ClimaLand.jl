# Canopy Autotrophic Respiration

This page documents the **autotrophic respiration** module used in ClimaLand's vegetation standalone runs. It follows the exact temperature‑response equation implemented in `src/standalone/Vegetation/autotrophic_respiration.jl`.

---

## Summary
- **Purpose:** compute autotrophic (plant) respiration of the canopy.
- **Inputs:** canopy (or leaf) temperature `T` and model parameters (`R_ref`, `Q10`, `T_ref`).
- **Core:** temperature scaling of a reference respiration rate using a fixed **Q10** formulation.
- **Outputs:** canopy autotrophic respiration rate `R_auto(T)` (e.g., W m⁻² or μmol CO₂ m⁻² s⁻¹ depending on configuration).

---

## Model formulation (as in the code)

Autotrophic respiration is calculated from a reference rate at a reference temperature using a **Q10** relationship:

```math
R_{\text{auto}}(T) \,=\, R_{\text{ref}}\; Q_{10}^{\;\dfrac{T - T_{\text{ref}}}{10}}\,.
```

Where:
- `R_ref` — respiration rate at the reference temperature `T_ref`,
- `Q10` — dimensionless temperature‑sensitivity parameter,
- `T` — canopy (or leaf) temperature (°C),
- `T_ref` — reference temperature (°C).

> Notes
> * Units of `R_auto` match the units of `R_ref`. The code maintains unit consistency; ensure your drivers and diagnostics use the same basis.
> * The module applies the formula directly; no additional acclimation or substrate‑dependence terms are included in this implementation.

---

## Variables & units

| Quantity | Symbol | Units | Notes |
|---|:--:|:--:|---|
| Temperature | `T` | °C | canopy or leaf temperature used by the module |
| Reference temperature | `T_ref` | °C | typically 20 or 25 |
| Reference respiration | `R_ref` | model‑dependent | rate at `T_ref` |
| Q10 parameter | `Q10` | – | temperature sensitivity |
| Autotrophic respiration | `R_auto(T)` | same as `R_ref` | output rate |

---

## Parameters

| Parameter | Description | Typical value/range |
|---|---|---|
| `T_ref` | reference temperature | 20–25 °C |
| `Q10` | temperature sensitivity | 1.5–2.5 |
| `R_ref` | base respiration at `T_ref` | calibrated/site‑specific |

---

