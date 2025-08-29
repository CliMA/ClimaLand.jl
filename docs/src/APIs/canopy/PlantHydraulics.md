# PlantHydraulics

```@meta
CurrentModule = ClimaLand.PlantHydraulics
```
## Models and Parameters

```@docs
ClimaLand.PlantHydraulics.PlantHydraulicsModel
ClimaLand.PlantHydraulics.PlantHydraulicsParameters
```

## Plant Hydraulics Parameterizations

```@docs
ClimaLand.PlantHydraulics.Weibull
ClimaLand.PlantHydraulics.LinearRetentionCurve
ClimaLand.PlantHydraulics.PrescribedSiteAreaIndex
ClimaLand.PlantHydraulics.AbstractTranspiration
ClimaLand.PlantHydraulics.DiagnosticTranspiration
```

## Constructor Methods

```@docs
ClimaLand.Canopy.PlantHydraulicsModel{FT}(
    domain,
    LAI::AbstractTimeVaryingInput,
    toml_dict::CP.AbstractTOMLDict;
    n_stem::Int = 0,
    n_leaf::Int = 1,
    h_stem::FT = FT(0),
    h_leaf::FT = FT(1),
    SAI::FT = toml_dict["SAI"],
    RAI::FT = toml_dict["RAI"],
    ai_parameterization = PlantHydraulics.PrescribedSiteAreaIndex{FT}(
        LAI,
        SAI,
        RAI,
    ),
    ν::FT = toml_dict["plant_nu"],
    S_s::FT = toml_dict["plant_S_s"], # m3/m3/MPa to m3/m3/m
    conductivity_model = PlantHydraulics.Weibull(toml_dict),
    retention_model = PlantHydraulics.LinearRetentionCurve(toml_dict),
    rooting_depth = clm_rooting_depth(domain.space.surface),
    transpiration = PlantHydraulics.DiagnosticTranspiration{FT}(),
) where {FT <: AbstractFloat}
```

## Plant Hydraulics Diagnostic Variables

```@docs
ClimaLand.PlantHydraulics.water_flux
ClimaLand.PlantHydraulics.effective_saturation
ClimaLand.PlantHydraulics.augmented_liquid_fraction
ClimaLand.PlantHydraulics.water_retention_curve
ClimaLand.PlantHydraulics.inverse_water_retention_curve
ClimaLand.PlantHydraulics.root_water_flux_per_ground_area!
ClimaLand.PlantHydraulics.hydraulic_conductivity
ClimaLand.PlantHydraulics.root_distribution
```
