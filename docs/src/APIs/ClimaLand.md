# ClimaLand

```@meta
CurrentModule = ClimaLand
```
## Integrated Land Model Types and methods

```@docs
ClimaLand.LandModel
ClimaLand.LandModel{FT}()
ClimaLand.LandModel{FT}(
    forcing,
    LAI,
    toml_dict::CP.AbstractTOMLDict,
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
    toml_dict::CP.AbstractTOMLDict,
    domain::Union{ClimaLand.Domains.Column, ClimaLand.Domains.SphericalShell};
) where {FT}
ClimaLand.LandHydrology
ClimaLand.LandHydrology{FT}()
ClimaLand.LandSoilBiogeochemistry
ClimaLand.LandSoilBiogeochemistry{FT}()
ClimaLand.SoilSnowModel
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
ClimaLand.PrognosticSoilConditions
```

## LandSoilBiogeochemistry

```@docs
ClimaLand.PrognosticMet
```
