# Photosynthesis

```@meta
CurrentModule = ClimaLand.Canopy
```

## Models and Parameters

```@docs
ClimaLand.Canopy.PModel
ClimaLand.Canopy.PModelParameters
ClimaLand.Canopy.PModelConstants
ClimaLand.Canopy.FarquharModel
ClimaLand.Canopy.FarquharParameters
ClimaLand.Canopy.OptimalityFarquharModel
ClimaLand.Canopy.OptimalityFarquharParameters
```

## FarquharModel Methods

```@docs
ClimaLand.Canopy.arrhenius_function
ClimaLand.Canopy.intercellular_co2
ClimaLand.Canopy.co2_compensation
ClimaLand.Canopy.rubisco_assimilation
ClimaLand.Canopy.light_assimilation
ClimaLand.Canopy.max_electron_transport
ClimaLand.Canopy.electron_transport
ClimaLand.Canopy.net_photosynthesis
ClimaLand.Canopy.optimality_max_photosynthetic_rates
ClimaLand.Canopy.moisture_stress
ClimaLand.Canopy.dark_respiration
ClimaLand.Canopy.compute_GPP
ClimaLand.Canopy.MM_Kc
ClimaLand.Canopy.MM_Ko
ClimaLand.Canopy.compute_Vcmax
```

## PModel Methods
ClimaLand.Canopy.compute_full_pmodel_outputs
ClimaLand.Canopy.set_historical_cache!
ClimaLand.Canopy.update_optimal_EMA
ClimaLand.Canopy.make_PModel_callback
ClimaLand.Canopy.intrinsic_quantum_yield
ClimaLand.Canopy.compute_viscosity_ratio
ClimaLand.Canopy.compute_Kmm
ClimaLand.Canopy.optimal_co2_ratio_c3
ClimaLand.Canopy.intercellular_co2_pmodel
ClimaLand.Canopy.gs_co2_pmodel
ClimaLand.Canopy.gs_h2o_pmodel
ClimaLand.Canopy.vcmax_pmodel
ClimaLand.Canopy.compute_LUE
ClimaLand.Canopy.compute_mj_with_jmax_limitation
ClimaLand.Canopy.electron_transport_pmodel
ClimaLand.Canopy.co2_compensation_p
ClimaLand.Canopy.quadratic_soil_moisture_stress
ClimaLand.Canopy.compute_APAR
