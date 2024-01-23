struct AnalyticTimeVaryingInput{F <: Function} <: AbstractTimeVaryingInput
    # func here as to be GPU-compatible (e.g., splines are not)
    func::F
end

function TimeVaryingInput(input::Function; method = nothing, device = nothing)
    isnothing(method) ||
        @warn "Interpolation method is ignored for analytical functions"
    return AnalyticTimeVaryingInput(input)
end

function evaluate!(dest, input::AnalyticTimeVaryingInput, time)
    dest .= input.func(time)
    return nothing
end
