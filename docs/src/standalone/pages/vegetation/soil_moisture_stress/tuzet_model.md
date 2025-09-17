# Tuzet Soil Moisture (Water Potential) Stress

This page documents the **Tuzet** stress function used to scale leaf net assimilation and stomatal conductance as a function of plant water potential. It provides a smooth, sigmoidal decline in conductance with decreasing (more negative) leaf/xylem water potential. See [Tuzet2003](@citet) and [Duursma2012](@citet). The ClimaLand code can be found [here](https://github.com/CliMA/ClimaLand.jl/blob/main/src/standalone/Vegetation/soil_moisture_stress.jl).

---

## Summary
- **Purpose:** compute a water‑status stress multiplier $\beta \in [0,1]$ from leaf (or xylem) water potential.
- **Inputs:** leaf/xylem water potential $\psi$ (Pa), parameters controlling the midpoint and slope.
- **Outputs:** scalar stress factor $\beta$ applied within stomatal/photosynthesis sub‑models.

---

## Model formulation

We compute the stress factor using the following formulation:

```math
\beta = \min\Bigg(1, \frac{1 + e^{s_c \, p_c}}{1 + e^{s_c \, (p_c - p_{leaf})}}\Bigg)
```

with

The leaf water pressure $p_{leaf}$ (Pa),

The reference water pressure $p_{c}$ (Pa),

The sensitivity to low water pressure $s_{c}$ (Pa$^{-1}$).

---

## Variables & units

| Quantity | Symbol | Units | Notes |
|---|:--:|:--:|---|
| Leaf water pressure | $p_{leaf}$ | Pa | typically ≤ 0 (tension) |
| Reference potential | $p_{c}$ | Pa | location of curve midpoint |
| Sensitivity | $s_{c}$ | Pa⁻¹ | slope/steepness of response |
| Stress multiplier | $\beta$ | – | 0 (full stress) … 1 (no stress) |

---
