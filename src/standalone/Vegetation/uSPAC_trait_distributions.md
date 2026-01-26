# uSPAC Stomatal Conductance with Trait Distributions

## Overview

This document describes the mathematical implementation of the uSPAC (unitless Soil Plant Atmosphere Continuum) stomatal conductance model with **climate-dependent trait distributions** in ClimaLand.jl.

**Key features:**
- Mechanistic stomatal conductance using Π-group dimensional analysis
- Trait heterogeneity via Gaussian quadrature integration over 3D trait space
- Climate-dependent variance: dry climates → high diversity (niche partitioning), wet → low diversity (competition)
- Captures emergent properties: compensatory dynamics, temporal asynchrony, non-linear buffering

---

## Part 1: Core uSPAC Model (Mean-Field)

### Dimensional Analysis Framework

The uSPAC model reduces the full hydraulic system to four dimensionless Π-groups that control stomatal conductance:

#### Π-groups Definition

```julia
ΠR = |ψ_g50| / |ψ_x50|      # Isohydry index (0-1)
ΠF = E₀ / (kx · ψ_g50)       # Flux control
ΠT = (K_SR · ψ_s) / E₀       # Supply capacity
ΠS = ψ_g50 / ψ_s             # Suitability
```

Where:
- `ψ_g50`: Leaf water potential at 50% stomatal closure (MPa)
- `ψ_x50` (P50): Xylem water potential at 50% conductance loss (MPa)
- `kx`: Xylem hydraulic conductance (m day⁻¹ MPa⁻¹)
- `E₀`: Reference evapotranspiration (m day⁻¹)
- `K_SR`: Soil-root hydraulic conductance
- `ψ_s`: Soil water potential (MPa)

#### From Π-groups to Stomatal Conductance

**Step 1: Compute fractional water-water point (fww)**
```julia
fww = 1 - (1/(2ΠR)) × [1 + ΠF/2 - √((1 + ΠF/2)² - 2ΠF·ΠR)]
```

**Step 2: Compute critical saturation thresholds**
```julia
# Wilting point (β = 0.05 × fww)
s_w = [(ΠT/(2β_w·ΠS)) × (√(1 + 4β_w·ΠS²·(2(1-β_w) - β_w·ΠF)/(ΠT·(1-(1-β_w)·ΠR))) - 1)]^(-1/b)

# Stomatal closure point (β = 0.95 × fww)
s* = [(ΠT/(2β*·ΠS)) × (√(1 + 4β*·ΠS²·(2(1-β*) - β*·ΠF)/(ΠT·(1-(1-β*)·ΠR))) - 1)]^(-1/b)
```

**Step 3: Water stress factor β(s)**
```julia
r = clamp((s - s_w) / (s* - s_w), 0, 1)
β = fww × r
```

**Step 4: Stomatal conductance**
```julia
E = β × E₀                           # Actual ET (m day⁻¹)
gsw_leaf = (ρw/Mw) × E × (P_air/VPD) # Leaf-level (mol m⁻² s⁻¹)
gsw_canopy = gsw_leaf × LAI          # Canopy-level (ground area basis)
```

### Climate-Dependent Mean Traits

Mean traits vary with aridity following calibrated β coefficients:

```julia
# Aridity normalization (0 = dry, 1 = wet)
aridity_norm = aridity_index / 1000  # P/ET₀ index

# P50: Cavitation resistance (always negative)
μ_P50 = -exp(βψx50_base - βψx50_slope × aridity_norm)
# Dry → -10 MPa (very resistant), Wet → -0.5 MPa (vulnerable)

# kx: Hydraulic conductance (coordinated with P50)
μ_logkx = βkx_base + βkx_coord × log(-μ_P50)
# Safety-efficiency trade-off (Liu et al. 2019)

# ΠR: Isohydry index (0-1)
μ_ΠR = 1 / (1 + exp(-(βΠR_base - βΠR_slope × aridity_norm)))
# Dry → 1 (anisohydric), Wet → 0 (isohydric)
```

---

## Part 2: Trait Distribution Extension

### Motivation: Why Heterogeneity Matters

**Jensen's Inequality:** For non-linear functions, E[f(x)] ≠ f(E[x])

The uSPAC model is highly non-linear:
- Exponential vulnerability curves (P50 thresholds)
- Piecewise β(s) response
- Π-group combinations involve products and square roots

**Mean-field assumption misses:**
1. **Compensatory dynamics:** Drought-sensitive plants close → water for resistant plants
2. **Temporal asynchrony:** Staggered P50 thresholds → gradual ecosystem transitions
3. **Non-linear buffering:** Heterogeneity prevents catastrophic collapse

### 3D Trait Space

Traits follow a multivariate normal distribution in transformed space:

```julia
p(log(kx), P50, ΠR) = MVN(μ, Σ)

μ = [μ_logkx, μ_P50, μ_ΠR]  # Climate-dependent means

Σ = [σ_logkx²              σ_logkx·σ_P50·ρ_kx_P50    σ_logkx·σ_ΠR·ρ_kx_ΠR  ]
    [σ_logkx·σ_P50·ρ_kx_P50  σ_P50²                   σ_P50·σ_ΠR·ρ_P50_ΠR  ]
    [σ_logkx·σ_ΠR·ρ_kx_ΠR    σ_P50·σ_ΠR·ρ_P50_ΠR       σ_ΠR²                ]
```

**Physical bounds enforced:**
- P50 < 0 (always negative water potential)
- 0 ≤ ΠR ≤ 1 (bounded isohydry index)
- kx > 0 (positive conductance)

### Gaussian-Hermite Quadrature Integration

Instead of sampling millions of trait combinations, we use **optimal quadrature nodes**:

```julia
E[gsw_leaf] = ∫∫∫ gsw_leaf(logkx, P50, ΠR) · p(logkx, P50, ΠR) d(logkx) dP50 dΠR

            ≈ Σᵢ wᵢ · gsw_leaf(logkxᵢ, P50ᵢ, ΠRᵢ)
```

**Quadrature properties:**
- n = 3 → 3³ = 27 points in 3D (default)
- n = 5 → 125 points (higher accuracy, slower)
- n = 7 → 343 points (research mode)
- Weights wᵢ sum to 1.0
- Exact for polynomials up to degree 2n-1

**Node placement:** Transformed from standard Gauss-Hermite tables via Cholesky decomposition of Σ = LL^T

---

## Part 3: Climate-Dependent Variance (Phase 1.5)

### Ecological Hypothesis

```
    Trait Variance (σ)          Trait Coordination (ρ)
    
2.4 |  🌵 ●──────────         0.7 |              ●─── 🌳
    |     │                       |            ╱
2.0 |     │                   0.6 |          ╱
    |     │                       |        ╱
1.6 |     └──────● 🌳          0.5 |  ●───╯ 🌵
    └──────────────────           └──────────────────
      Dry        Wet                Dry        Wet
```

**Dry climates (aridity_norm < 0.3):**
- HIGH variance (σ_P50 ≈ 2.4 MPa)
- WEAK coordination (ρ_kx_P50 ≈ 0.5)
- Mechanism: **Niche partitioning** for limited water
  - Deep-rooted + resistant xylem (P50 = -8 MPa)
  - Shallow-rooted + opportunistic (P50 = -1 MPa)
  - Diverse water-use strategies coexist
  
**Wet climates (aridity_norm > 0.7):**
- LOW variance (σ_P50 ≈ 1.6 MPa)
- STRONG coordination (ρ_kx_P50 ≈ 0.7)
- Mechanism: **Competitive exclusion**
  - All plants converge to optimal trait combination
  - Strong safety-efficiency trade-off
  - Homogeneous canopy

### Mathematical Implementation

```julia
# Compute aridity stress: 0 (wet) → 1 (dry)
aridity_stress = 1 - aridity_norm

# Standard deviations increase with aridity stress
σ_logkx = σ_logkx_base + α_σ_logkx × aridity_stress
σ_P50 = σ_P50_base + α_σ_P50 × aridity_stress
σ_ΠR = σ_ΠR_base + α_σ_ΠR × aridity_stress

# Correlations weaken in dry climates
ρ_kx_P50 = ρ_kx_P50_base + α_ρ_kx_P50 × aridity_stress
```

**Default parameters (best-guess from literature):**
```julia
# Base values (wet climate reference)
σ_logkx_base = 0.5     # Liu et al. 2019
σ_P50_base = 1.5       # Choat et al. 2012
σ_ΠR_base = 0.15
ρ_kx_P50_base = 0.7    # Manzoni et al. 2013

# Climate-dependent modifiers
α_σ_logkx = 0.3        # +60% variance in driest climates
α_σ_P50 = 1.0          # +67% variance in driest climates
α_σ_ΠR = 0.1
α_ρ_kx_P50 = -0.2      # -29% coordination in driest climates
```

**Gradient across climates:**
| Aridity Index | Climate | σ_P50 (MPa) | ρ_kx_P50 |
|--------------|---------|-------------|----------|
| 0.0 (dry) | Desert/Savanna | 2.5 | 0.50 |
| 0.2 | Semi-arid | 2.3 | 0.54 |
| 0.5 | Mediterranean | 2.0 | 0.60 |
| 0.7 | Temperate | 1.8 | 0.64 |
| 0.9 (wet) | Rainforest | 1.6 | 0.68 |

---

## Part 4: Emergent Ecosystem Properties

### Drought Response Comparison

**Scenario:** Soil water potential drops from -0.5 MPa to -6 MPa over 30 days

#### 🌳 Wet Climate Ecosystem (Low σ_P50 = 1.6 MPa)

All plants have similar P50 ≈ -2 MPa

```
Day  0: ████████████████████  100% active
Day 10: ████████████████████  100% (still above threshold)
Day 15: ███████              35% (SHARP DROP - most hit P50)
Day 20: ██                   10% (ecosystem collapse)
Day 30: ▌                    5%  (near-total shutdown)
```

**Characteristics:**
- Sharp, catastrophic transition
- Low buffering capacity
- Synchronized stomatal closure
- High vulnerability to drought

#### 🌵 Dry Climate Ecosystem (High σ_P50 = 2.4 MPa)

Diverse traits: P50 ranges from -1 to -8 MPa

```
Day  0: ████████████████████  100% active
Day 10: ████████████████▌     75% (opportunists P50=-1 close)
Day 15: ████████████          60% (moderates P50=-3 close)
Day 20: ████████              40% (intermediates P50=-5 struggle)
Day 25: █████                 25% (only resistant left)
Day 30: ███                   15% (deep-rooted P50=-8 persist)
```

**Characteristics:**
- Smooth, gradual transition
- High buffering capacity
- Asynchronous stomatal closure
- Resilience through diversity

### Compensatory Dynamics

When drought-sensitive plants close stomata:
1. Reduced transpiration → slower soil drying
2. More water available for resistant plants
3. Resistant plants sustain ecosystem function longer
4. **Result:** Total canopy conductance declines smoothly, not catastrophically

---

## Part 5: Implementation in ClimaLand

### Parameter Struct

```julia
Base.@kwdef struct uSPACPiParameters{FT <: AbstractFloat}
    # Mean trait calibration (aridity-dependent)
    βkx_base::FT      # Log-scale intercept for kx
    βkx_coord::FT     # Coordination with P50
    βψx50_base::FT    # P50 exponential base
    βψx50_slope::FT   # P50 aridity slope
    βΠR_base::FT      # ΠR logistic base
    βΠR_slope::FT     # ΠR aridity slope
    
    # Trait variance (climate-dependent)
    σ_logkx_base::FT = FT(0.5)
    σ_P50_base::FT = FT(1.5)
    σ_ΠR_base::FT = FT(0.15)
    α_σ_logkx::FT = FT(0.3)
    α_σ_P50::FT = FT(1.0)
    α_σ_ΠR::FT = FT(0.1)
    
    # Trait correlations (climate-dependent)
    ρ_kx_P50_base::FT = FT(0.7)
    ρ_kx_ΠR_base::FT = FT(-0.3)
    ρ_P50_ΠR_base::FT = FT(-0.2)
    α_ρ_kx_P50::FT = FT(-0.2)
    
    # Quadrature settings
    n_quad::Int = 3                      # 3, 5, or 7
    use_trait_distribution::Bool = false  # Enable heterogeneity
    
    # Physical parameters
    b::FT = FT(4.38)        # Soil moisture exponent
    β_star_frac::FT = FT(0.95)
    β_w_frac::FT = FT(0.05)
    gsw_max::FT = FT(Inf)
end
```

### Update Loop (Simplified)

```julia
function update_canopy_conductance!(p, Y, model::uSPACConductancePi, canopy)
    pars = model.parameters
    
    # Extract state: soil saturation s, E₀, VPD, etc.
    aridity_norm = extract_aridity(...)
    
    # Compute mean traits from climate
    μ_P50 = -exp(pars.βψx50_base - pars.βψx50_slope × aridity_norm)
    μ_logkx = pars.βkx_base + pars.βkx_coord × log(-μ_P50)
    μ_ΠR = 1 / (1 + exp(-(pars.βΠR_base - pars.βΠR_slope × aridity_norm)))
    
    if pars.use_trait_distribution
        # Climate-dependent variance
        aridity_stress = 1 - aridity_norm
        σ_P50 = pars.σ_P50_base + pars.α_σ_P50 × aridity_stress
        σ_logkx = pars.σ_logkx_base + pars.α_σ_logkx × aridity_stress
        ρ_kx_P50 = pars.ρ_kx_P50_base + pars.α_ρ_kx_P50 × aridity_stress
        
        # Generate quadrature
        quad = generate_trait_quadrature(
            FT, pars.n_quad,
            μ_logkx, μ_P50, μ_ΠR,
            σ_logkx, σ_P50, σ_ΠR,
            ρ_kx_P50, ρ_kx_ΠR, ρ_P50_ΠR
        )
        
        # Integrate over traits
        gsw_leaf = zero(s)
        for i in 1:length(quad.samples)
            logkx_i, P50_i, ΠR_i = quad.samples[i]
            kx_i = exp(logkx_i)
            
            # Compute Π-groups for this trait sample
            ΠF_i, ΠT_i, ΠS_i = pi_groups_from_traits(kx_i, P50_i, ΠR_i, ...)
            
            # uSPAC algebra: Π → (fww, s*, sw) → β(s) → gsw
            fww_i = uspac_fww(ΠR_i, ΠF_i)
            sstar_i = uspac_s_of_beta(0.95 × fww_i, ...)
            sw_i = uspac_s_of_beta(0.05 × fww_i, ...)
            β_i = fww_i × clamp((s - sw_i)/(sstar_i - sw_i), 0, 1)
            gsw_i = (ρw/Mw) × β_i × E₀ × (P_air / VPD)
            
            # Accumulate weighted contribution
            gsw_leaf += quad.weights[i] × gsw_i
        end
    else
        # Mean-field (original, faster)
        ΠF, ΠT, ΠS = pi_groups_from_traits(μ_kx, μ_P50, μ_ΠR, ...)
        fww = uspac_fww(μ_ΠR, ΠF)
        # ... standard computation
        gsw_leaf = (ρw/Mw) × β × E₀ × (P_air / VPD)
    end
    
    # Convert to canopy resistance
    r_stomata_canopy = 1 / (gsw_leaf × LAI × R × T_air / P_air)
end
```

---

## Part 6: Usage Examples

### Mean-Field Mode (Default, Fast)

```julia
params = uSPACPiParameters(
    βkx_base = -6.5,
    βkx_coord = -2.5,
    βψx50_base = 2.3,
    βψx50_slope = 0.8,
    βΠR_base = -1.5,
    βΠR_slope = 2.0,
    use_trait_distribution = false  # Mean traits only
)

# Same speed as before, no trait integration
```

### Trait Distribution Mode (Research, ~27× Slower)

```julia
params = uSPACPiParameters(
    βkx_base = -6.5,
    βkx_coord = -2.5,
    βψx50_base = 2.3,
    βψx50_slope = 0.8,
    βΠR_base = -1.5,
    βΠR_slope = 2.0,
    # Climate-dependent variance
    σ_P50_base = 1.5,      # Wet climate (Choat et al. 2012)
    α_σ_P50 = 1.0,         # Dry climates: σ → 2.5 MPa
    ρ_kx_P50_base = 0.7,   # Wet climate (Manzoni et al. 2013)
    α_ρ_kx_P50 = -0.2,     # Dry climates: ρ → 0.5
    n_quad = 3,            # 27 points in 3D
    use_trait_distribution = true  # ENABLE
)

# Integrates over 27 trait samples per gridcell per timestep
# Captures emergent properties (compensatory dynamics, buffering)
```

---

## Part 7: Testable Predictions

### 1. Remote Sensing Heterogeneity
**Prediction:** Spatial variance increases with aridity stress

```julia
σ_NDVI(aridity) ∝ σ_P50(aridity) = σ_base + α × aridity_stress
```

**Test:** Landsat 30m NDVI variance vs. MAP/PET gradient
- Dry ecosystems: High sub-grid heterogeneity
- Wet ecosystems: Homogeneous canopy

### 2. FLUXNET Temporal Variance
**Prediction:** Flux variance during drought correlates with aridity

```julia
CV_GPP(aridity) = σ_GPP / μ_GPP ∝ σ_P50(aridity)
```

**Test:** Compare coefficient of variation across AmeriFlux sites
- Dry sites (CA-NS, US-Wkg): High σ_GPP during summer
- Wet sites (US-Ha1, BR-Sa3): Low σ_GPP, sharp transitions

### 3. TRY Database Validation
**Prediction:** Within-site trait variance follows aridity gradient

**Test:** Extract P50, kx measurements grouped by site
- Compute σ_P50 for each site
- Correlate with site aridity index
- Expected: σ_P50 ∝ (1 - aridity_norm)

### 4. Ecosystem Resilience
**Prediction:** Dry ecosystems buffer drought stress better

**Metrics:**
- GPP decline rate: Dry → gradual, Wet → sharp
- Recovery time: Dry → fast (resistant plants maintained), Wet → slow (full regrowth)

**Test:** MODIS GPP trajectories during 2012 US drought, 2018 European heatwave
- Classify sites by aridity
- Compare time to 50% GPP reduction
- Expected: Dry sites take 2-3× longer to reach threshold

---

## Part 8: Performance Characteristics

| Mode | Speed | Memory | Accuracy | Use Case |
|------|-------|--------|----------|----------|
| Mean-field | 1× | 1× | Mean traits only | Production runs, benchmarking |
| Trait distribution (n=3) | 27× | 27× | Captures heterogeneity | Research, variance studies |
| Trait distribution (n=5) | 125× | 125× | Higher accuracy | Sensitivity analysis |
| Trait distribution (n=7) | 343× | 343× | Research only | Method validation |

**Optimization opportunities (future):**
- GPU acceleration (trait samples are embarrassingly parallel)
- Adaptive quadrature (skip integration in low-variance regions)
- Cached Π-group computations across similar trait samples

---

## Part 9: Key References

**Dimensional analysis & uSPAC:**
- Bassiouni et al. (2023): Universal stomatal behavior from Π-groups
- Sperry et al. (2017): Predicting stomatal responses via hydraulics

**Trait coordination:**
- Liu et al. (2019): kx-P50 safety-efficiency trade-off
- Manzoni et al. (2013): Coordination shapes drought response (ρ_kx_P50 ≈ 0.7)

**Trait variance:**
- Choat et al. (2012): Global P50 database (σ ≈ 1.5 MPa across species)
- Anderegg et al. (2016): Trait diversity buffers drought mortality

**Quadrature methods:**
- Golub & Welsch (1969): Gauss-Hermite quadrature tables
- Xiu & Karniadakis (2002): Polynomial chaos for uncertainty quantification
