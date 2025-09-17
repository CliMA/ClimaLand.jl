# Piecewise Soil Moisture Stress

This page documents the **piecewise** (linear/threshold with curvature) soil moisture stress function used in ClimaLand's vegetation module. The stress factor $\beta \in [0,1]$ scales leaf photosynthesis (and thus stomatal conductance) as soil moisture declines, following a simple, interpretable formulation. See [Egea2011](@citet), the ClimaLand code is [here](https://github.com/CliMA/ClimaLand.jl/blob/main/src/standalone/Vegetation/soil_moisture_stress.jl).

---

## Summary
- **Purpose:** convert soil wetness to a *stress multiplier* that reduces assimilation and stomatal conductance.
- **Inputs:** volumetric soil water content $\theta$ (or relative water content), model thresholds.
- **Outputs:** scalar stress factor $\beta$ applied within the photosynthesis/stomatal sub-models.

---

## Model formulation

Let $\theta$ be the volumetric water content at some depth. Define two thresholds:

First, $\theta_{low}$ — **wilting** (or residual) water content: below this, stress is total ($\beta = 0$).

Second, $\theta_{high}$ — **Field capacity (or porosity)** water content: above this, no stress ($\beta = 1$).

With $\theta_{low} < \theta_{high}$, the stress factor is

```math
\beta(\theta) =
\begin{cases}
0, & \theta \leq \theta_{low}, \\
\left( \dfrac{\theta - \theta_{low}}{\theta_{high} - \theta_{low}} \right)^c, & \theta_{low} < \theta < \theta_{high}, \\
1, & \theta \geq \theta_{high}. \\
\end{cases}
```

---

## Variables & units

| Quantity | Symbol | Units | Notes |
|---|:--:|:--:|---|
| Volumetric water content | $\theta$ | m³ m⁻³ | root-zone or layer-weighted |
| Low water content threshold | $\theta_{low}$ | m³ m⁻³ | wilting point or residual  |
| High water content threshold | $\theta_{high}$ | m³ m⁻³ | field capacity or porosity |
| Stress multiplier | $\beta$ | – | 0 (full stress) … 1 (no stress) |
| Curvature parameter | $c$ | – | controls concavity of the stress response |

---

## Parameters

| Parameter | Description | Example range |
|---|---|---|
| $\theta_{low}$ | lower threshold (wilting or residual) | 0.05–0.15 |
| $\theta_{high}$ | upper threshold (field capacity or porosity) | 0.55–0.75 |
| $c$ | curvature parameter shaping the transition | 1–5 |

---

## Numerical details & coupling

Soil moisture, $\theta$, is either prescribed as a single value in the root zone (for standalone canopy models, with a prescribed soil component), or else we compute $\beta$ as a factor of depth using $\theta$ as a function of depth, as specified by the prognostic soil model. In the latter case, we then average $\beta(z)$ over the column using the root distribution function (also a function of$z$), as a weighting factor.

Scalar stress, $\beta$, multiplies leaf assimilation, $A_{n}$, which in turn reduces $g_{s}$ that scales with $A_{n}$.

---