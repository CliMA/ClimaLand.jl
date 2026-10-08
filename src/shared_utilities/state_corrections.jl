export LogLinearFactor

"""
    LogLinearFactor{FT, N}

A bounded multiplicative correction factor of `N` state variables `x`,

    f(x) = exp(clamp(b + Σᵢ wᵢ (xᵢ - μᵢ) / σᵢ, -L, L)),

with bias `b`, weights `w`, feature centers `μ` and scales `σ`, and the bound
`L = log_fmax` on the logarithm of the factor. With `b = 0` and `w = 0` the
factor is one.

Such factors carry empirical corrections, for example regressions of the
residuals of a model against flux-tower observations, into the
parameterizations: they scale a conductance, an albedo, or a clumping index
inside the existing flux solves, so the corrected model conserves energy and
water exactly as the uncorrected one does, and the correction cannot leave the
range `[1/f_max, f_max]` where the state leaves the range of the training data.

The factor is called with the feature values as positional arguments,
`f(x₁, …, x_N)`, and can be broadcast over fields.
$(DocStringExtensions.FIELDS)
"""
struct LogLinearFactor{FT <: AbstractFloat, N}
    "Bias b of the logarithm of the factor"
    bias::FT
    "Weights w of the standardized features"
    weights::NTuple{N, FT}
    "Feature centers μ"
    center::NTuple{N, FT}
    "Feature scales σ (> 0)"
    scale::NTuple{N, FT}
    "Bound L on |log f|"
    log_fmax::FT
end

"""
    LogLinearFactor{FT}(;
        weights,
        center = zeros(length(weights)),
        scale = ones(length(weights)),
        bias = 0,
        log_fmax = log(3),
    )

Construct a [`LogLinearFactor`](@ref) from vectors or tuples of weights,
centers and scales.
"""
function LogLinearFactor{FT}(;
    weights,
    center = ntuple(_ -> zero(FT), length(weights)),
    scale = ntuple(_ -> one(FT), length(weights)),
    bias = zero(FT),
    log_fmax = log(FT(3)),
) where {FT}
    N = length(weights)
    @assert length(center) == N && length(scale) == N
    @assert all(>(0), scale)
    return LogLinearFactor{FT, N}(
        FT(bias),
        NTuple{N, FT}(weights),
        NTuple{N, FT}(center),
        NTuple{N, FT}(scale),
        FT(log_fmax),
    )
end

@inline function (f::LogLinearFactor{FT, N})(x::Vararg{Any, N}) where {FT, N}
    z = f.bias
    z += sum(
        ntuple(i -> f.weights[i] * (x[i] - f.center[i]) / f.scale[i], Val(N)),
    )
    return exp(clamp(z, -f.log_fmax, f.log_fmax))
end

Base.Broadcast.broadcastable(f::LogLinearFactor) = Ref(f)
Base.eltype(::LogLinearFactor{FT}) where {FT} = FT
