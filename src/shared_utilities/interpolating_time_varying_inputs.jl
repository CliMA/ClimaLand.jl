InterpolatingTimeVaryingInput =
    Union{InterpolatingTimeVaryingInput0D, InterpolatingTimeVaryingInput2D}

"""
    NearestNeighbor

Return the value corresponding to the point closest to the input time.
"""
struct NearestNeighbor <: AbstractInterpolationMethod end

"""
    LinearInterpolation

Perform linear interpolation between the two neighboring points.
"""
struct LinearInterpolation <: AbstractInterpolationMethod end

"""
    in(time, itp::InterpolatingTimeVaryingInput0D)

Check if the given `time` is in the range of definition for `itp`.
"""
function Base.in(time, itp::InterpolatingTimeVaryingInput0D)
    return itp.range[1] <= time <= itp.range[2]
end
