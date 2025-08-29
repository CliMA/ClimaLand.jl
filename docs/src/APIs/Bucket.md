# Bucket

```@meta
CurrentModule = ClimaLand.Bucket
```
## Types and Constructors

```@docs
ClimaLand.Bucket.BucketModelParameters
ClimaLand.Bucket.BucketModelParameters(
    ::Type{FT};
) where {FT <: AbstractFloat}
ClimaLand.Bucket.BucketModelParameters(
    toml_dict::CP.AbstractTOMLDict;
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
