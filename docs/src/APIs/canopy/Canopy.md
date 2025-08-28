# Canopy

```@meta
CurrentModule = ClimaLand.Canopy
```
## Canopy Model and Parameters

```@docs
ClimaLand.Canopy.CanopyModel
ClimaLand.Canopy.CanopyModel{FT}(;
    autotrophic_respiration::AbstractAutotrophicRespirationModel{FT},
    radiative_transfer::AbstractRadiationModel{FT},
    photosynthesis::AbstractPhotosynthesisModel{FT},
    conductance::AbstractStomatalConductanceModel{FT},
    hydraulics::AbstractPlantHydraulicsModel{FT},
    energy = PrescribedCanopyTempModel{FT}(),
    sif = Lee2015SIFModel{FT}(),
    boundary_conditions::B,
    parameters::SharedCanopyParameters{FT, PSE},
    domain::Union{
        ClimaLand.Domains.Point,
        ClimaLand.Domains.Plane,
        ClimaLand.Domains.SphericalSurface,
    },
) where {FT, B, PSE}
ClimaLand.Canopy.CanopyModel{FT}(
    domain::Union{
        ClimaLand.Domains.Point,
        ClimaLand.Domains.Plane,
        ClimaLand.Domains.SphericalSurface,
    },
    forcing::NamedTuple,
    LAI::AbstractTimeVaryingInput,
    toml_dict::CP.AbstractTOMLDict;
    z_0m = toml_dict["canopy_momentum_roughness_length"],
    z_0b = toml_dict["canopy_scalar_roughness_length"],
    prognostic_land_components = (:canopy,),
    autotrophic_respiration = AutotrophicRespirationModel{FT}(),
    radiative_transfer = TwoStreamModel{FT}(domain),
    photosynthesis = FarquharModel{FT}(domain),
    conductance = MedlynConductanceModel{FT}(domain),
    hydraulics = PlantHydraulicsModel{FT}(domain, LAI, toml_dict),
    energy = BigLeafEnergyModel{FT}(),
    sif = Lee2015SIFModel{FT}(),
) where {FT}
ClimaLand.Canopy.SharedCanopyParameters
ClimaLand.Canopy.AbstractCanopyComponent
```

## Canopy Model Boundary Fluxes

```@docs
ClimaLand.Canopy.AbstractCanopyBC
ClimaLand.Canopy.AtmosDrivenCanopyBC
```
