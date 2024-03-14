InterpolatingTimeVaryingInput =
    Union{InterpolatingTimeVaryingInput0D, InterpolatingTimeVaryingInput2D}

"""
    in(time, itp::InterpolatingTimeVaryingInput)

Check if the given `time` is in the range of definition for `itp`.
"""
function Base.in(time, itp::InterpolatingTimeVaryingInput)
    return itp.range[1] <= time <= itp.range[2]
end
