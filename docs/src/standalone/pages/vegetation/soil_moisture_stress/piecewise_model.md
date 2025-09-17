# Piecewise Soil Moisture Stress

This page documents the **piecewise** (linear/threshold) soil moisture stress function used in ClimaLand's vegetation module. The stress factor `S_θ ∈ [0,1]` scales photosynthetic and/or stomatal processes as soil moisture declines, following a simple, interpretable piecewise formulation.

---

## Summary
- **Purpose:** convert soil wetness to a *stress multiplier* that reduces assimilation, stomatal conductance, and/or mesophyll conductance.
- **Inputs:** volumetric soil water content `θ` (or relative water content), model thresholds.
- **Outputs:** scalar stress factor `S_θ` applied within the photosynthesis/stomatal sub‑models.

---

## Model formulation

Let `θ` be the volumetric water content of the plant‑available root zone (or a weighted layer average if roots span multiple layers). Define two thresholds:

- `θ_w` — **wilting** (or residual) water content: below this, stress is total (`S_θ = 0`).
- `θ_c` — **critical** water content: above this, no stress (`S_θ = 1`).

With `θ_w < θ_c`, the stress factor is

```math
S_\theta(\theta) =
\begin{cases}
0, & \theta \le \theta_w, \\
\dfrac{\theta - \theta_w}{\theta_c - \theta_w}, & \theta_w < \theta < \theta_c, \\
1, & \theta \ge \theta_c. \\
\end{cases}
```

Optionally, a **minimum floor** `S_{min}` can be enforced so that `S_θ ≥ S_{min}` for numerical stability.

---

## Variables & units

| Quantity | Symbol | Units | Notes |
|---|:--:|:--:|---|
| Volumetric water content | `θ` | m³ m⁻³ | root‑zone or layer‑weighted |
| Wilting water content | `θ_w` | m³ m⁻³ | often ≈ residual θ |
| Critical water content | `θ_c` | m³ m⁻³ | often near field capacity or a calibrated point |
| Stress multiplier | `S_θ` | – | 0 (full stress) … 1 (no stress) |

---

## Parameters

| Parameter | Description | Typical value/range |
|---|---|---|
| `θ_w` | lower threshold (wilting) | 0.05–0.15 |
| `θ_c` | upper threshold (critical) | 0.15–0.35 |
| `Smin` | minimum stress floor | 0.0–0.1 (optional) |
| `root_weights` | weights for multi‑layer averaging | normalized over layers |

---

## Numerical details & coupling
- `θ` is typically computed from hydrology as a (possibly weighted) root‑zone average.
- `S_θ` multiplies stomatal/biochemical capacities (e.g., `gₛ`, `V_{cmax}`, `J_{max}`) depending on configuration.
- Smoothing at the breakpoints can be enabled (e.g., short cubic ramp) if required for differentiability.

---

