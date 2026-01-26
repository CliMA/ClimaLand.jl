"""
    minimum_finite([f=identity], A; kwargs...)

Calculate `minimum(f, A)` while ignoring any values that are not finite, e.g., `Inf` or
`NaN`.

If `A` is a colorant array with multiple channels (e.g., `Array{RGB}`), the `min` comparison
is done in channel-wise sense.

The supported `kwargs` are those of `minimum(f, A; kwargs...)`.
"""
function minimum_finite(f, A::AbstractArray{T}; kwargs...) where T
    # FIXME(johnnychen94): if `typeof(f(first(A))) != eltype(A)`, this function is not type-stable.
    mapreduce(IfElse(isfinite, f, typemax), minc, A; kwargs...)
end
minimum_finite(A::AbstractArray; kwargs...) = minimum_finite(identity, A; kwargs...)

"""
    maximum_finite([f=identity], A; kwargs...)

Calculate `maximum(f, A)` while ignoring any values that are not finite, e.g., `Inf` or
`NaN`.

If `A` is a colorant array with multiple channels (e.g., `Array{RGB}`), the `max` comparison
is done in channel-wise sense.

The supported `kwargs` are those of `maximum(f, A; kwargs...)`.
"""
function maximum_finite(f, A::AbstractArray{T}; kwargs...) where T
    # FIXME(johnnychen94): if `typeof(f(first(A))) != eltype(A)`, this function is not type-stable
    mapreduce(IfElse(isfinite, f, typemin), maxc, A; kwargs...)
end
maximum_finite(A::AbstractArray; kwargs...) = maximum_finite(identity, A; kwargs...)

"""
    sumfinite([f=identity], A; kwargs...)

Compute `sum(f, A)` while ignoring any non-finite values.

The supported `kwargs` are those of `sum(f, A; kwargs...)`.
"""
sumfinite(A; kwargs...) = sumfinite(identity, A; kwargs...)

if Base.VERSION >= v"1.1"
    sumfinite(f, A; kwargs...) = sum(IfElse(isfinite, f, zero), A; kwargs...)
else
    sumfinite(f, A; kwargs...) = sum(IfElse(isfinite, f, zero).(A); kwargs...)
end

"""
    meanfinite([f=identity], A; kwargs...)

Compute `mean(f, A)` while ignoring any non-finite values.

The supported `kwargs` are those of `sum(f, A; kwargs...)`.
"""
meanfinite(A; kwargs...) = meanfinite(identity, A; kwargs...)

if Base.VERSION >= v"1.1"
    function meanfinite(f, A; kwargs...)
        s = sumfinite(f, A; kwargs...)
        n = sum(IfElse(isfinite, x->true, x->false), A; kwargs...)   # TODO: replace with `Returns`
        return s./n
    end
else
    function meanfinite(f, A; kwargs...)
        s = sumfinite(f, A; kwargs...)
        n = sum(IfElse(isfinite, x->true, x->false).(A); kwargs...)
        return s./n
    end
end

"""
    varfinite(A; kwargs...)

Compute the variance of `A`, ignoring any non-finite values.

The supported `kwargs` are those of `sum(f, A; kwargs...)`.

!!! note
    This function can produce a seemingly suprising result if the input array is an RGB
    image. To make it more clear, the implementation is made so that
    `varfinite(img) ≈ varfinite(RGB.(img))` holds for any gray-scale image. See also
    https://github.com/JuliaGraphics/ColorVectorSpace.jl#abs-and-abs2 for more information.

"""
function varfinite end

if Base.VERSION >= v"1.1"
    function varfinite(A; kwargs...)
        m = meanfinite(A; kwargs...)
        n = sum(IfElse(isfinite, x->true, x->false), A; kwargs...)   # TODO: replace with `Returns`
        s = sum(IfElse(isfinite, identity, zero), abs2.(A .- m); kwargs...)
        return s ./ max.(0, (n .- 1))
    end
else
    function varfinite(A; kwargs...)
        m = meanfinite(A; kwargs...)
        n = sum(IfElse(isfinite, x->true, x->false).(A); kwargs...)
        s = sum(IfElse(isfinite, identity, zero).(abs2.(A .- m)); kwargs...)
        return s ./ max.(0, (n .- 1))
    end
end
