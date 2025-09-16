# ClimaLand

```@meta
CurrentModule = ClimaLand
```
## Integrated Land Model Types and methods

```@docs
ClimaLand.LandModel
ClimaLand.LandModel{FT}(
    forcing,
    LAI,
    toml_dict::CP.ParamDict,
    domain::Union{
        ClimaLand.Domains.Column,
        ClimaLand.Domains.SphericalShell,
        ClimaLand.Domains.HybridBox,
    },
    Δt;
) where {FT}
ClimaLand.SoilCanopyModel
ClimaLand.SoilCanopyModel{FT}(
    forcing,
    LAI,
    toml_dict::CP.ParamDict,
    domain::Union{ClimaLand.Domains.Column, ClimaLand.Domains.SphericalShell};
) where {FT}
ClimaLand.LandSoilBiogeochemistry
ClimaLand.LandSoilBiogeochemistry{FT}(
    forcing,
    toml_dict::CP.ParamDict,
    domain::Union{ClimaLand.Domains.Column, ClimaLand.Domains.SphericalShell};
) where {FT}
ClimaLand.SoilSnowModel
ClimaLand.SoilSnowModel{FT}(
    forcing,
    toml_dict::CP.ParamDict,
    domain::Union{
        ClimaLand.Domains.Column,
        ClimaLand.Domains.SphericalShell,
        ClimaLand.Domains.HybridBox,
    },
    Δt;
) where {FT}
ClimaLand.land_components
ClimaLand.lsm_aux_vars
ClimaLand.lsm_aux_types
ClimaLand.lsm_aux_domain_names
```

## Land Hydrology

```@docs
ClimaLand.infiltration_capacity
ClimaLand.infiltration_at_point
ClimaLand.PrognosticRunoff
ClimaLand.RunoffBC
```

## SoilCanopyModel

```@docs
ClimaLand.RootExtraction
```

## LandSoilBiogeochemistry

```@docs
ClimaLand.PrognosticMet
```
