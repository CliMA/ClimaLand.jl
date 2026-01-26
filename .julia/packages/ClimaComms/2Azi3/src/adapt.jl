import Adapt

"""
    Adapt.adapt_structure(::Type{<:AbstractArray}, context::AbstractCommsContext)

Adapt a given context to a context with a device associated with the given array type.

# Example

```julia
Adapt.adapt_structure(Array, ClimaComms.context(ClimaComms.CUDADevice())) -> ClimaComms.CPUSingleThreaded()
```

!!! note
    By default, adapting to `Array` creates a `CPUSingleThreaded` device, and
    there is currently no way to conver to a CPUMultiThreaded device.
"""
Adapt.adapt_structure(to::Type{<:AbstractArray}, ctx::AbstractCommsContext) =
    context(Adapt.adapt(to, device(ctx)))

"""
    Adapt.adapt_structure(::Type{<:AbstractArray}, device::AbstractDevice)

Adapt a given device to a device associated with the given array type.

# Example

```julia
Adapt.adapt_structure(Array, ClimaComms.CUDADevice()) -> ClimaComms.CPUSingleThreaded()
```

!!! note
    By default, adapting to `Array` creates a `CPUSingleThreaded` device, and
    there is currently no way to conver to a CPUMultiThreaded device.
"""
Adapt.adapt_structure(::Type{<:AbstractArray}, device::AbstractDevice) =
    CPUSingleThreaded()
