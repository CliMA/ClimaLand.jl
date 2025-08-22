# ClimaLand

```@meta
CurrentModule = ClimaLand
```
## Integrated Land Model Types and methods

```@docs
ClimaLand.LandModel
ClimaLand.LandModel{FT}(;
    soilco2_type::Type{MM},
    soilco2_args::NamedTuple = (;),
    land_args::NamedTuple = (;),
    soil_model_type::Type{SM},
    soil_args::NamedTuple = (;),
    canopy_component_types::NamedTuple = (;),
    canopy_component_args::NamedTuple = (;),
    canopy_model_args::NamedTuple = (;),
    snow_model_type::Type{SnM},
    snow_args::NamedTuple = (;),
) where {
    FT,
    SM <: Soil.EnergyHydrology{FT},
    MM <: Soil.Biogeochemistry.SoilCO2Model{FT},
    SnM <: Snow.SnowModel,
}
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
    soil = Soil.EnergyHydrology{FT}(
        domain,
        forcing,
        toml_dict;
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
        additional_sources = (ClimaLand.RootExtraction{FT}(),),
    ),
    soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(
        domain,
        Soil.Biogeochemistry.SoilDrivers(
            Soil.Biogeochemistry.PrognosticMet(soil.parameters),
            PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5)),
            forcing.atmos,
        ),
    ),
    canopy = Canopy.CanopyModel{FT}(
        Domains.obtain_surface_domain(domain),
        (;
            atmos = forcing.atmos,
            radiation = forcing.radiation,
            ground = ClimaLand.PrognosticGroundConditions{FT}(),
        ),
        LAI,
        toml_dict;
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
    ),
    snow = Snow.SnowModel(
        FT,
        ClimaLand.Domains.obtain_surface_domain(domain),
        forcing,
        toml_dict,
        Δt;
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
    ),
) where {FT}
ClimaLand.SoilCanopyModel
ClimaLand.SoilCanopyModel{FT}(
    forcing,
    LAI,
    toml_dict::CP.AbstractTOMLDict,
    domain::Union{ClimaLand.Domains.Column, ClimaLand.Domains.SphericalShell};
    soil = Soil.EnergyHydrology{FT}(
        domain,
        forcing,
        toml_dict;
        prognostic_land_components = (:canopy, :soil, :soilco2),
        additional_sources = (ClimaLand.RootExtraction{FT}(),),
    ),
    soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(
        domain,
        Soil.Biogeochemistry.SoilDrivers(
            Soil.Biogeochemistry.PrognosticMet(soil.parameters),
            PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5)),
            forcing.atmos,
        ),
    ),
    canopy = Canopy.CanopyModel{FT}(
        Domains.obtain_surface_domain(domain),
        (;
            atmos = forcing.atmos,
            radiation = forcing.radiation,
            ground = ClimaLand.PrognosticSoilConditions{FT}(),
        ),
        LAI,
        toml_dict;
        prognostic_land_components = (:canopy, :soil, :soilco2),
    ),
) where {FT}
ClimaLand.LandHydrology
ClimaLand.LandHydrology{FT}(;
    land_args::NamedTuple = (;),
    soil_model_type::Type{SM},
    soil_args::NamedTuple = (;),
    surface_water_model_type::Type{SW},
    surface_water_args::NamedTuple = (;),
) where {
    FT,
    SM <: Soil.AbstractSoilModel{FT},
    SW <: Pond.AbstractSurfaceWaterModel{FT},
}
ClimaLand.LandSoilBiogeochemistry
ClimaLand.LandSoilBiogeochemistry{FT}(;
    land_args::NamedTuple,
    soil_args::NamedTuple = (;),
    soilco2_args::NamedTuple = (;),
) where {FT}
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
