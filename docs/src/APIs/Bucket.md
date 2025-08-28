# Bucket

```@meta
CurrentModule = ClimaLand.Bucket
```
## Types and Constructors

```@docs
ClimaLand.Bucket.BucketModelParameters
ClimaLand.Bucket.BucketModelParameters(
    ::Type{FT};
    albedo,
    z_0m,
    z_0b,
    τc,
    kwargs...,
) where {FT <: AbstractFloat}
ClimaLand.Bucket.BucketModelParameters(
    toml_dict::CP.AbstractTOMLDict;
    albedo,
    W_f = toml_dict["land_bucket_capacity"],
    f_bucket = toml_dict["bucket_capacity_fraction"],
    z_0b = toml_dict["bucket_z_0b"],
    ρc_soil = toml_dict["bucket_soil_heat_capacity"],
    f_snow = toml_dict["critical_snow_fraction"],
    z_0m = toml_dict["bucket_z_0m"],
    σS_c = toml_dict["critical_snow_water_equivalent"],
    p = toml_dict["bucket_beta_decay_exponent"],
    κ_soil = toml_dict["bucket_soil_conductivity"],
    τc = toml_dict["tau_c"],
)
ClimaLand.Bucket.PrescribedBaregroundAlbedo
ClimaLand.Bucket.PrescribedSurfaceAlbedo
ClimaLand.Bucket.BucketModel
```

## Misc Functions

```@docs
ClimaLand.Bucket.surface_albedo
ClimaLand.Bucket.beta_factor
```

