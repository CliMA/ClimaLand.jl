# Tuzet Soil Moisture (Water Potential) Stress

This page documents the **Tuzet** stress function used to scale stomatal (and/or photosynthetic) capacity as a function of plant water potential. It provides a smooth, sigmoidal decline in conductance with decreasing (more negative) leaf/xylem water potential.

---

## Summary
- **Purpose:** compute a water‑status stress multiplier `S_ψ ∈ [0,1]` from leaf (or xylem) water potential.
- **Inputs:** leaf/xylem water potential `ψ` (MPa), parameters controlling the midpoint and slope.
- **Outputs:** scalar stress factor `S_ψ` applied within stomatal/photosynthesis sub‑models.

---

## Model formulation

Following Tuzet et al. (2003), the water‑potential response is

```math
S_\psi(\psi) = \frac{1}{1 + \exp\big[a\,(\psi - \psi_{50})\big]},
```

where:

- `ψ` — leaf or xylem water potential (MPa; negative in tension),
- `ψ_{50}` — midpoint potential where `S_ψ = 0.5`,
- `a` — dimensionless slope (steepness). Larger `a` yields a sharper transition.

Some implementations include a residual floor `S_{min}` so that `S_ψ = S_{min} + (1-S_{min})/(1 + exp[a(ψ-ψ_{50})])`.

> **Choice of `ψ`:** Depending on configuration, `ψ` may be the soil water potential (root‑weighted), stem/xylem, or leaf potential derived from a hydraulic model.

---

## Variables & units

| Quantity | Symbol | Units | Notes |
|---|:--:|:--:|---|
| Water potential | `ψ` | MPa | typically ≤ 0 |
| Midpoint potential | `ψ_{50}` | MPa | where `S_ψ = 0.5` |
| Slope | `a` | – | logistic steepness |
| Stress multiplier | `S_ψ` | – | 0 (full stress) … 1 (no stress) |

---

## Parameters

| Parameter | Description | Typical value/range |
|---|---|---|
| `ψ50` | midpoint water potential | −2 to −1 MPa (species/site dependent) |
| `a` | slope parameter | 2–12 |
| `Smin` | residual floor | 0.0–0.1 (optional) |
| `state_variable` | which potential to use (`:leaf`, `:xylem`, `:soil`) | `:leaf` (common) |

---

## Numerical details & coupling
- `ψ` is read from the hydraulics/soil module at each step (or derived diagnostically).
- `S_ψ` multiplies `gₛ` and, optionally, biochemical capacities.
- The logistic form is smooth/differentiable; no special smoothing is required.

---
