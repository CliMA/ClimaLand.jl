# [Reference](@id function_reference)

## Spatial transformation

- `restrict` originally lives in [ImageTransformation.jl](https://github.com/JuliaImages/ImageTransformations.jl). This function is still reexported by ImageTransformation.

```@docs
restrict
```

## Discrete gradient operator

```@autodocs
Modules = [ImageBase.FiniteDiff]
order = [:module, :function]
```

## Statistics

```@docs
minimum_finite
maximum_finite
meanfinite
varfinite
sumfinite
```
