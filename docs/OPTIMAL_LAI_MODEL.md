# Optimal LAI Model Implementation

## Overview

This document describes the implementation of the prognostic leaf area index (LAI) model based on Zhou et al. (2025), "A General Model for the Seasonal to Decadal Dynamics of Leaf Area," published in Global Change Biology.

## Model Description

The optimal LAI model dynamically predicts LAI by balancing photosynthetic carbon assimilation against the cost of maintaining leaf area. The model updates LAI at local noon using an exponential moving average (EMA) to account for acclimation timescales.

### Key Equations

1. **Top-of-canopy assimilation rate** (inverting Beer-Lambert law):
   ```
   Ao = A / (1 - exp(-k*L))
   ```
   where `A` is canopy-integrated assimilation, `k` is extinction coefficient, `L` is LAI.

2. **Maximum LAI constraint** (ensures bottom leaves are viable):
   ```
   L_max = -log(z/(k*Ao)) / k
   ```
   where `z` is the minimum assimilation threshold.

3. **Optimal LAI** (balances benefit vs. cost):
   Solves: `L/μ - 1 + exp(-k*L) = 0`
   where `μ = m*Ao` combines maintenance cost `m` and assimilation capacity `Ao`.

4. **EMA update** (accounts for acclimation):
   ```
   L_new = (1 - α) * min(L_opt, L_max) + α * L_current
   ```
   where `α` is the memory parameter (higher = more weight on history).

## Parameters

Three parameters control the model behavior (defined in `toml/default_parameters.toml`):

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `optimal_lai_m` | 0.3 | unitless | Maintenance cost parameter. Represents marginal cost of additional leaf area. |
| `optimal_lai_z` | 12.227 | μmol CO₂ m⁻² s⁻¹ | Minimum assimilation threshold. Limits LAI in low-light conditions. |
| `optimal_lai_alpha` | 0.933 | unitless | Acclimation timescale. α = 1 - 1 day/τ where τ ≈ 15 days. |

## Implementation Details

### Current Architecture

The optimal LAI model is integrated into the P-model photosynthesis:

- **Location**: `src/standalone/Vegetation/pmodel.jl`
- **Auxiliary Variable**: `L` (leaf area index) stored in `p.canopy.photosynthesis.L`
- **Update Frequency**: Daily at local noon
- **Update Function**: `update_optimal_LAI(L, A, cosθs, canopy, local_noon_mask)`
- **Callback**: Integrated into `make_PModel_callback`

### Core Functions

All LAI computation functions are in `pmodel.jl`:

```julia
compute_Ao(A, k, L)           # Top-of-canopy assimilation
compute_L_max(k, z, Ao)       # Maximum allowable LAI
g(μ, k, L)                    # Optimality condition
dgdL(μ, k, L)                 # Derivative for Newton-Raphson
compute_L_opt(μ, k, L)        # Optimal LAI via Newton-Raphson
compute_L(L, Ao, m, k, α, z, mask)  # EMA update at local noon
```

### OptimalLAIModel Structure

For future standalone use, `src/standalone/Vegetation/optimal_lai.jl` provides:

```julia
OptimalLAIParameters{FT}      # Parameter structure
OptimalLAIModel{FT}            # Model structure (for future use)
```

This structure is currently not actively used but is available for future integration as a separate canopy component.

## Usage

### In Simulations

When using `PModel` for photosynthesis, LAI is automatically predicted:

```julia
# Set up photosynthesis with P-model
photosynthesis = PModel{FT}(surface_domain, toml_dict; is_c3)

# The canopy model will automatically update LAI at local noon
canopy = Canopy.CanopyModel{FT}(
    surface_domain,
    forcing,
    LAI,  # This can be initial LAI or prescribed LAI
    toml_dict;
    photosynthesis,
    # ... other components
)
```

### Diagnostics

The predicted LAI is available as diagnostic variable:

```julia
output_vars = ["lai", "lai_pred", "gpp", ...]
```

- `"lai"`: Prescribed LAI (from biomass model or observations)
- `"lai_pred"`: Predicted LAI (from optimal LAI model)

### Example: Ozark FLUXNET Site

See `experiments/integrated/fluxnet/ozark_pmodel.jl` for a complete example:

```julia
# Model setup automatically includes optimal LAI
photosynthesis = PModel{FT}(surface_domain, toml_dict; is_c3)

# Request both prescribed and predicted LAI in diagnostics
output_vars = ["gpp", "shf", "lhf", "swu", "lwu", "swc", "swe", "tsoil", 
               "lai", "lai_pred"]
```

## Customizing Parameters

To modify the default parameters, update `toml/default_parameters.toml`:

```toml
["optimal_lai_m"]
value = 0.25  # Lower cost → higher LAI

["optimal_lai_z"]
value = 10.0  # Lower threshold → can grow more LAI in shade

["optimal_lai_alpha"]
value = 0.8   # Lower α → faster acclimation (τ ≈ 5 days)
```

## Initialization

LAI is initialized in `set_historical_cache!`:

```julia
L = p.canopy.photosynthesis.L
L .= 1.0  # Initial LAI of 1.0 m²/m²
```

The model then adjusts LAI based on environmental conditions over the first acclimation period (~15 days by default).

## References

Zhou, B., Cai, W., Zhu, Z., Wang, H., Harrison, S. P., & Prentice, C. (2025). A General Model for the Seasonal to Decadal Dynamics of Leaf Area. *Global Change Biology*, 31(1), e70125. https://doi.org/10.1111/gcb.70125

## Code Structure

```
src/standalone/Vegetation/
├── pmodel.jl              # Core LAI functions + PModel integration
├── optimal_lai.jl         # OptimalLAIModel structure (future use)
└── Canopy.jl              # Includes both files

toml/
└── default_parameters.toml  # Default parameter values

experiments/integrated/fluxnet/
└── ozark_pmodel.jl        # Example usage

src/diagnostics/
├── define_diagnostics.jl  # Diagnostic variable definitions
└── land_compute_methods.jl  # Compute method: lai_pred → p.canopy.photosynthesis.L
```

## Future Extensions

The `OptimalLAIModel` structure in `optimal_lai.jl` is designed to support:

1. **Standalone LAI component**: Could be added as separate field in `CanopyModel`
2. **Multiple LAI models**: Switch between prescribed, optimal, or other LAI models
3. **Custom parameterizations**: Use `OptimalLAIParameters` for spatially-varying parameters
4. **Separate callbacks**: Independent LAI update timing from photosynthesis

These extensions would require modifications to the `CanopyModel` structure but the core functionality is already implemented.
