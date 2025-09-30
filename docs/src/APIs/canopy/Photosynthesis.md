# Photosynthesis

```@meta
CurrentModule = ClimaLand.Canopy
```

## Models and Parameters

```@docs
ClimaLand.Canopy.PModel
ClimaLand.Canopy.PModel(domain, toml_dict::CP.ParamDict)
ClimaLand.Canopy.PModelParameters
ClimaLand.Canopy.PModelConstants
ClimaLand.Canopy.FarquharModel
ClimaLand.Canopy.FarquharModel{FT}(
    domain,
    toml_dict::CP.ParamDict;
) where {FT <: AbstractFloat}
ClimaLand.Canopy.FarquharParameters
ClimaLand.Canopy.FarquharParameters(toml_dict::CP.ParamDict)
```

## FarquharModel Methods

```@docs
ClimaLand.Canopy.arrhenius_function
ClimaLand.Canopy.intercellular_co2_farquhar
ClimaLand.Canopy.co2_compensation_farquhar
ClimaLand.Canopy.rubisco_assimilation_farquhar
ClimaLand.Canopy.light_assimilation_farquhar
ClimaLand.Canopy.max_electron_transport_farquhar
ClimaLand.Canopy.electron_transport_farquhar
ClimaLand.Canopy.net_photosynthesis
ClimaLand.Canopy.dark_respiration_farquhar
ClimaLand.Canopy.MM_Kc
ClimaLand.Canopy.MM_Ko
ClimaLand.Canopy.compute_Vcmax_farquhar
```

## PModel Methods

```@docs
ClimaLand.Canopy.compute_full_pmodel_outputs
ClimaLand.Canopy.set_historical_cache!
ClimaLand.Canopy.update_optimal_EMA
ClimaLand.Canopy.make_PModel_callback
ClimaLand.Canopy.intrinsic_quantum_yield
ClimaLand.Canopy.compute_viscosity_ratio
ClimaLand.Canopy.compute_Kmm
ClimaLand.Canopy.intercellular_co2_pmodel
ClimaLand.Canopy.gs_co2_pmodel
ClimaLand.Canopy.gs_h2o_pmodel
ClimaLand.Canopy.compute_mj
ClimaLand.Canopy.compute_mc
ClimaLand.Canopy.vcmax_pmodel
ClimaLand.Canopy.compute_LUE
ClimaLand.Canopy.compute_mj_with_jmax_limitation
ClimaLand.Canopy.electron_transport_pmodel
ClimaLand.Canopy.co2_compensation_pmodel
ClimaLand.Canopy.quadratic_soil_moisture_stress
ClimaLand.Canopy.compute_APAR_canopy_moles
ClimaLand.Canopy.compute_APAR_leaf_moles
```

## Shared Methods

```@docs
ClimaLand.Canopy.gross_photosynthesis
```
