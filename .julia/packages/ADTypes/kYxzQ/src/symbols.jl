"""
    ADTypes.Auto(package::Symbol)
    ADTypes.Auto(nothing)::NoAutoDiff

A shortcut that converts an AD package name into an instance of [`AbstractADType`](@ref), with all parameters set to their default values.

!!! warning

    This function is type-unstable by design and might lead to suboptimal performance.
    In most cases, you should never need it: use the individual backend types directly.

# Example

```jldoctest
import ADTypes
backend = ADTypes.Auto(:Zygote)

# output

ADTypes.AutoZygote()
```
"""
Auto(package::Symbol, args...; kws...) = Auto(Val(package), args...; kws...)

for backend in (:ChainRules, :Diffractor, :Enzyme, :FastDifferentiation,
    :FiniteDiff, :FiniteDifferences, :ForwardDiff, :Mooncake, :PolyesterForwardDiff,
    :ReverseDiff, :Symbolics, :Tapir, :Tracker, :Zygote)
    @eval Auto(::Val{$(QuoteNode(backend))}, args...; kws...) = $(Symbol(:Auto, backend))(
        args...; kws...)
end

Auto(::Nothing) = NoAutoDiff()
