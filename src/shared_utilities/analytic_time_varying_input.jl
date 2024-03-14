struct AnalyticTimeVaryingInput{F <: Function} <: AbstractTimeVaryingInput
    # func here has to be GPU-compatible (e.g., splines are not) and reasonably fast (e.g.,
    # no large allocations)
    func::F
end

# _kwargs... is needed to seamlessly support the other TimeVaryingInputs.
function TimeVaryingInput(input::Function; method = nothing, _kwargs...)
    isnothing(method) ||
        @warn "Interpolation method is ignored for analytical functions"
    return AnalyticTimeVaryingInput(input)
end

function evaluate!(
    dest,
    input::AnalyticTimeVaryingInput,
    time,
    args...;
    kwargs...,
)
    dest .= input.func(time, args...; kwargs...)
    return nothing
end
