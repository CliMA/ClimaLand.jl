"""
    InterpolatingTimeVaryingInput0D

The constructor for InterpolatingTimeVaryingInput0D is not supposed to be used directly, unless you
know what you are doing. The constructor does not perform any check and does not take care of
GPU compatibility. It is responsibility of the user-facing constructor TimeVaryingInput() to do so.

`times` and `vales` may have different float types, but they must be the same length, and we
assume that they have been sorted to be monotonically increasing in time, without repeated
values for the same timestamp.
"""
struct InterpolatingTimeVaryingInput0D{
    AA1 <: AbstractArray,
    AA2 <: AbstractArray,
    M <: AbstractInterpolationMethod,
    CC <: ClimaComms.AbstractCommsContext,
    R <: Tuple,
} <: AbstractTimeVaryingInput
    # AA1 and AA2 could be different because of different FTs

    """Independent coordinate"""
    times::AA1

    """Variable"""
    vals::AA2

    """Interpolation method"""
    method::M

    """ClimaComms context"""
    context::CC

    """Range of times over which the interpolator is defined. range is always defined on the
    CPU. Used by the in() function."""
    range::R
end

# GPU compatibility
function Adapt.adapt_structure(to, itp::InterpolatingTimeVaryingInput0D)
    times = Adapt.adapt_structure(to, itp.times)
    vals = Adapt.adapt_structure(to, itp.vals)
    method = Adapt.adapt_structure(to, itp.method)
    range = Adapt.adapt_structure(to, itp.range)
    # On a GPU, we have a "ClimaCore.DeviceSideContext"
    InterpolatingTimeVaryingInput0D(
        times,
        vals,
        method,
        DeviceSideContext(),
        range,
    )
end

function evaluate!(
    destination,
    itp::InterpolatingTimeVaryingInput0D,
    time,
    args...;
    kwargs...,
)
    time in itp || error("TimeVaryingInput does not cover time $time")
    if ClimaComms.device(itp.context) isa ClimaComms.CUDADevice
        CUDA.@cuda evaluate!(parent(destination), itp, time, itp.method)
    else
        evaluate!(parent(destination), itp, time, itp.method)
    end
    return nothing
end

function TimeVaryingInput(
    times::AbstractArray,
    vals::AbstractArray;
    method = LinearInterpolation(),
    context = ClimaComms.context(),
)
    issorted(times) || error("Can only interpolate with sorted times")
    length(times) == length(vals) ||
        error("times and vals have different lengths")

    # When device is CUDADevice, ArrayType will be a CUDADevice, so that times and vals get
    # copied to the GPU.
    ArrayType = ClimaComms.array_type(ClimaComms.device(context))

    range = (times[begin], times[end])
    return InterpolatingTimeVaryingInput0D(
        ArrayType(times),
        ArrayType(vals),
        method,
        context,
        range,
    )
end

function evaluate!(
    dest,
    itp::InterpolatingTimeVaryingInput0D,
    time,
    ::NearestNeighbor,
)
    # Nearest neighbor interpolation: just pick the values corresponding to the entry in
    # itp.times that is closest to the given time.

    index = searchsortednearest(itp.times, time)

    dest .= itp.vals[index]

    return nothing
end

"""
    evaluate!(
        dest,
        itp::InterpolatingTimeVaryingInput0D,
        time,
        ::LinearInterpolation,
        )

Write to `dest` the result of a linear interpolation of `itp` on the given `time`.
"""
function evaluate!(
    dest,
    itp::InterpolatingTimeVaryingInput0D,
    time,
    ::LinearInterpolation,
)
    # We perform linear interpolation with the two values that bracket the given time.
    # itp.times is sorted, so we find the first index that is after `time`, the previous
    # index is the other side of the bracket. searchsortedfirst returns the index of the
    # first element that is >= than the given. Also, given that we check that range[1] <=
    # time <= range[2], index will always be 1 <= index <= length(itp.times), so we have to
    # worry about the edge case where time == itp.times (because it returns 1). In that
    # case, we just return the value of vals[1] (we are on a node, no need for
    # interpolation).

    indep_vars = itp.times
    indep_value = time
    dep_vars = itp.vals
    dest .= linear_interpolation(indep_vars, dep_vars, indep_value)
    return nothing
end
