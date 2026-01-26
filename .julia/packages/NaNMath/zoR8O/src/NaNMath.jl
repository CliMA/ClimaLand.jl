module NaNMath

using OpenLibm_jll
const libm = OpenLibm_jll.libopenlibm

for f in (:sin, :cos, :tan, :asin, :acos, :acosh, :atanh,
          :log, :log2, :log10, :log1p, :lgamma)
    @eval begin
        function ($f)(x::Real)
            xf = float(x)
            x === xf && throw(MethodError($f, (x,)))
            ($f)(xf)
        end
        $(f !== :lgamma ? :(($f)(x) = (Base.$f)(x)) : :())
    end
end

Base.@assume_effects :total lgamma(x::Float64) = ccall(("lgamma",libm), Float64, (Float64,), x)
Base.@assume_effects :total lgamma(x::Float32) = ccall(("lgammaf",libm), Float32, (Float32,), x)

for f in (:sin, :cos, :tan)
    @eval begin
        function ($f)(x::T) where T<:Union{Float16, Float32, Float64}
            isinf(x) ? T(NaN) : (Base.$f)(x)
        end
    end
end

for f in (:asin, :acos, :atanh)
    @eval begin
        function ($f)(x::T) where T<:Union{Float16, Float32, Float64}
            abs(x) > T(1) ? T(NaN) : (Base.$f)(x)
        end
    end
end
function acosh(x::T) where T<:Union{Float16, Float32, Float64}
    x < T(1) ? T(NaN) : Base.acosh(x)
end

for f in (:log, :log2, :log10)
    @eval begin
        function ($f)(x::T) where T<:Union{Float16, Float32, Float64}
            x < 0 ? T(NaN) : (Base.$f)(x)
        end
    end
end

function log1p(x::T) where T<:Union{Float16, Float32, Float64}
    x < T(-1) ? T(NaN) : Base.log1p(x)
end

for f in (:sqrt,)
    @eval ($f)(x) = (Base.$f)(x)
end

for f in (:max, :min)
    @eval ($f)(x, y) = (Base.$f)(x, y)
end

sqrt(x::T) where {T<:Union{Float16, Float32, Float64}} = x < T(0) ? T(NaN) : Base.Intrinsics.sqrt_llvm(x)
sqrt(x::T) where {T<:AbstractFloat} = x < T(0) ? T(NaN) : Base.sqrt(x)
sqrt(x::Real) = sqrt(float(x))

# Don't override built-in ^ operator
Base.@assume_effects :total pow(x::Float64, y::Float64) = ccall((:pow,libm),  Float64, (Float64,Float64), x, y)
Base.@assume_effects :total pow(x::Float32, y::Float32) = ccall((:powf,libm), Float32, (Float32,Float32), x, y)
# We `promote` first before converting to floating pointing numbers to ensure that
# e.g. `pow(::Float32, ::Int)` ends up calling `pow(::Float32, ::Float32)`
pow(x::Real, y::Real) = pow(promote(x, y)...)
pow(x::T, y::T) where {T<:Real} = pow(float(x), float(y))
pow(x, y) = ^(x, y)

# The following combinations are safe, so we can fall back to ^
pow(x::Number, y::Integer) = x^y
pow(x::Real, y::Integer) = x^y
pow(x::Complex, y::Complex) = x^y

"""
NaNMath.sum(A)

##### Args:
* `A`: An array of floating point numbers

##### Returns:
*    Returns the sum of all elements in the array, ignoring NaN's.

##### Examples:
```julia
using NaNMath
NaNMath.sum([1., 2., NaN]) # result: 3.0
```
"""
function sum(x::AbstractArray{T}) where T<:AbstractFloat
    if length(x) == 0
        result = zero(eltype(x))
    else
        result = convert(eltype(x), NaN)
        for i in x
            if !isnan(i)
                if isnan(result)
                    result = i
                else
                    result += i
                end
            end
        end
    end

    if isnan(result)
        @warn "All elements of the array, passed to \"sum\" are NaN!"
    end
    return result
end

"""
NaNMath.median(A)

##### Args:
* `A`: An array of floating point numbers

##### Returns:
*   Returns the median of all elements in the array, ignoring NaN's.
    Returns NaN for an empty array or array containing NaNs only.

##### Examples:
```jldoctest
julia> using NaNMath

julia> NaNMath.median([1., 2., 3., NaN])
2.

julia> NaNMath.median([1., 2., NaN])
1.5

julia> NaNMath.median([NaN])
NaN
```
"""
median(x::AbstractArray{<:AbstractFloat}) = median(collect(Iterators.flatten(x)))

function median(x::AbstractVector{<:AbstractFloat})

    x = sort(filter(!isnan, x))

    n = length(x)
    if n == 0
        return convert(eltype(x), NaN)
    elseif isodd(n)
        ind = ceil(Int, n/2)
        return x[ind]
    else
        ind = Int(n/2)
        lower = x[ind]
        upper = x[ind+1]
        return (lower + upper) / 2
    end

end

"""
NaNMath.maximum(A)

##### Args:
* `A`: An array of floating point numbers

##### Returns:
*    Returns the maximum of all elements in the array, ignoring NaN's.

##### Examples:
```julia
using NaNMath
NaNMath.maximum([1., 2., NaN]) # result: 2.0
```
"""
function maximum(x::AbstractArray{T}) where T<:AbstractFloat
    result = convert(eltype(x), NaN)
    for i in x
        if !isnan(i)
            if (isnan(result) || i > result)
                result = i
            end
        end
    end
    return result
end

"""
NaNMath.minimum(A)

##### Args:
* `A`: An array of floating point numbers

##### Returns:
*    Returns the minimum of all elements in the array, ignoring NaN's.

##### Examples:
```julia
using NaNMath
NaNMath.minimum([1., 2., NaN]) # result: 1.0
```
"""
function minimum(x::AbstractArray{T}) where T<:AbstractFloat
    result = convert(eltype(x), NaN)
    for i in x
        if !isnan(i)
            if (isnan(result) || i < result)
                result = i
            end
        end
    end
    return result
end

"""
NaNMath.extrema(A)

##### Args:
* `A`: An array of floating point numbers

##### Returns:
*    Returns the minimum and maximum of all elements in the array, ignoring NaN's.

##### Examples:
```julia
using NaNMath
NaNMath.extrema([1., 2., NaN]) # result: 1.0, 2.0
```
"""
function extrema(x::AbstractArray{T}) where T<:AbstractFloat
    resultmin, resultmax = convert(eltype(x), NaN), convert(eltype(x), NaN)
    for i in x
        if !isnan(i)
            if (isnan(resultmin) || i < resultmin)
                resultmin = i
            end
            if (isnan(resultmax) || i > resultmax)
                resultmax = i
            end
        end
    end
    return resultmin, resultmax
end

"""
NaNMath.mean(A)

##### Args:
* `A`: An array of floating point numbers

##### Returns:
*    Returns the arithmetic mean of all elements in the array, ignoring NaN's.

##### Examples:
```julia
using NaNMath
NaNMath.mean([1., 2., NaN]) # result: 1.5
```
"""
function mean(x::AbstractArray{T}) where T<:AbstractFloat
    return mean_count(x)[1]
end

"""
Returns a tuple of the arithmetic mean of all elements in the array, ignoring NaN's,
and the number of non-NaN values in the array.
"""
function mean_count(x::AbstractArray{T}) where T<:AbstractFloat
    z = zero(eltype(x))
    sum = z
    count = 0
    @simd for i in x
        count += ifelse(isnan(i), 0, 1)
        sum += ifelse(isnan(i), z, i)
    end
    result = sum / count
    return (result, count)
end

"""
NaNMath.var(A)

##### Args:
* `A`: A one dimensional array of floating point numbers

##### Returns:
* Returns the sample variance of a vector A. The algorithm will return
  an estimator of the  generative distribution's variance under the
  assumption that each entry of v is an IID drawn from that generative
  distribution. This computation is  equivalent to calculating \\
  sum((v - mean(v)).^2) / (length(v) - 1). NaN values are ignored.

##### Examples:
```julia
using NaNMath
NaNMath.var([1., 2., NaN]) # result: 0.5
```
"""
function var(x::Vector{T}) where T<:AbstractFloat
    mean_val, n = mean_count(x)
    if !isnan(mean_val)
        sum_square = zero(eltype(x))
        for i in x
            if !isnan(i)
                sum_square += (i - mean_val)^2
            end
        end
        return sum_square / (n - one(eltype(x)))
    else
        return mean_val # NaN or NaN32
    end
end

"""
NaNMath.std(A)

##### Args:
* `A`: A one dimensional array of floating point numbers

##### Returns:
* Returns the standard deviation of a vector A. The algorithm will return
  an estimator of the  generative distribution's standard deviation under the
  assumption that each entry of v is an IID drawn from that generative
  distribution. This computation is  equivalent to calculating \\
  sqrt(sum((v - mean(v)).^2) / (length(v) - 1)). NaN values are ignored.

##### Examples:
```julia
using NaNMath
NaNMath.std([1., 2., NaN]) # result: 0.7071067811865476
```
"""
function std(x::Vector{T}) where T<:AbstractFloat
    return sqrt(var(x))
end

"""
    NaNMath.min(x, y)

Compute the IEEE 754-2008 compliant minimum of `x` and `y`. As of version 0.6 of Julia,
`Base.min(x, y)` will return `NaN` if `x` or `y` is `NaN`. `NanMath.min` favors values over
`NaN`, and will return whichever `x` or `y` is not `NaN` in that case.

## Examples

```julia
julia> NanMath.min(NaN, 0.0)
0.0

julia> NaNMath.min(1, 2)
1
```
"""
min(x::T, y::T) where {T<:AbstractFloat} = ifelse((y < x) | (signbit(y) > signbit(x)),
                                           ifelse(isnan(y), x, y),
                                           ifelse(isnan(x), y, x))

"""
    NaNMath.max(x, y)

Compute the IEEE 754-2008 compliant maximum of `x` and `y`. As of version 0.6 of Julia,
`Base.max(x, y)` will return `NaN` if `x` or `y` is `NaN`. `NaNMath.max` favors values over
`NaN`, and will return whichever `x` or `y` is not `NaN` in that case.

## Examples

```julia
julia> NaNMath.max(NaN, 0.0)
0.0

julia> NaNMath.max(1, 2)
2
```
"""
max(x::T, y::T) where {T<:AbstractFloat} = ifelse((y > x) | (signbit(y) < signbit(x)),
                                           ifelse(isnan(y), x, y),
                                           ifelse(isnan(x), y, x))

min(x::Real, y::Real) = min(promote(x, y)...)
max(x::Real, y::Real) = max(promote(x, y)...)

function min(x::BigFloat, y::BigFloat)
    isnan(x) && return y
    isnan(y) && return x
    return Base.min(x, y)
end

function max(x::BigFloat, y::BigFloat)
    isnan(x) && return y
    isnan(y) && return x
    return Base.max(x, y)
end

# Integers can't represent NaN
min(x::Integer, y::Integer) = Base.min(x, y)
max(x::Integer, y::Integer) = Base.max(x, y)

min(x::Real) = x
max(x::Real) = x

# Multi-arg versions
for f in (:min, :max)
    @eval ($f)(a, b, c, xs...) = Base.afoldl($f, ($f)(($f)(a, b), c), xs...)
end

# The functions `findmin`, `findmax`, `argmin`, and `argmax` are supported 
# to work correctly for the following iterable types:
_valtype(x::AbstractArray{<:AbstractFloat}) = eltype(x)
_valtype(x::Tuple{Vararg{AbstractFloat}}) = eltype(x)
_valtype(x::NamedTuple{<:Any, <:Tuple{Vararg{AbstractFloat}}}) = eltype(x)
_valtype(x::AbstractDict{<:Any,<:AbstractFloat}) = valtype(x)
_valtype(x) = error(
    "Iterables with value type AbstractFloat or its subtypes are supported.
    The provided input type $(typeof(x)) is not.
    Consider using the convert function before passing the iterable argument."

)

function _find_extreme(f,compare_op::Function, x)
    result_index = 1 # Note: default index value.
    result_value = convert(_valtype(x), NaN)

    for (k, v) in pairs(x)
        if !isnan(v)
            if (isnan(result_value) || compare_op(f(v),f(result_value)))
                result_index = k
                result_value = v
            end
        end
    end
    return f(result_value), result_index
end

"""
    NaNMath.findmin(f, domain) -> (f(x), index)

    NaNMath.findmin(domain) -> (x, index)

##### Args:
* `f`: A function applied to the elements of `domain`; 
  defaults to `identity` when `domain` is the only argument.
* `domain`: A non-empty collection of floating point numbers such that
  `f` is defined on elements of `domain`.

##### Returns:
* Returns a `Tuple` consisting of a value `f(x)` and the index of `x`
  in `domain`, ignoring NaN's, such that `f(x)` is minimized.
  If there are multiple minimal elements, then the first one will be returned.

If `domain` is a `NamedTuple` or dictionary-like `AbstractDict` L,
the function is applied to its values.  The returned index is a key `k`,
such that `f(L[k])` is minimized.

##### Examples:
```julia
julia> NaNMath.findmin([1., 1., 2., 2., NaN])
(1.0, 1)

julia> NaNMath.findmin(-, [1., 1., 2., 2., NaN])
(-2.0, 3)

julia> NaNMath.findmin(abs, Dict(:x => 3.0, :w => -2.2, :y => -3.0, :z => NaN))
(2.2, :w)
```
"""
function findmin end
findmin(f,x) = _find_extreme(f,<,x)
findmin(x) = findmin(identity,x)

"""
    NaNMath.findmax(f, domain) -> (f(x), index)

    NaNMath.findmax(domain) -> (x, index)

##### Args:
* `f`: A function applied to the elements of `domain`; 
  defaults to `identity` when `domain` is the only argument.
* `domain`: A non-empty collection of floating point numbers such that
  `f` is defined on elements of `domain`.

##### Returns:
* Returns a `Tuple` consisting of a value `f(x)` and the index of `x`
  in `domain`, ignoring NaN's, such that `f(x)` is maximized.
  If there are multiple maximal elements, then the first one will be returned.

If `domain` is a `NamedTuple` or dictionary-like `AbstractDict` L,
the function `f` is applied to its values.  The returned index is a key `k`,
such that `f(L[k])` is maximized.

##### Examples:
```julia
julia> NaNMath.findmax([1., 1., 2., 2., NaN])
(2.0, 3)

julia> NaNMath.findmax(-, [1., 1., 2., 2., NaN])
(-1.0, 1)

julia> NaNMath.findmax(abs, Dict(:x => 3.0, :w => -2.2, :y => -3.0, :z => NaN))
(3.0, :y)
```
"""
function findmax end
findmax(f,x) = _find_extreme(f,>,x)
findmax(x) = findmax(identity,x) 

"""
    NaNMath.argmin(f, domain) -> x

##### Args:
* `f`: A function applied to the elements of `domain`; 
  defaults to `identity` when `domain` is the only argument.
* `domain`: A non-empty collection of floating point numbers such that
  `f` is defined on elements of `domain`.

##### Returns:
* Returns a value `x` in `domain`, ignoring NaN's, for which `f(x)` is minimized.
  If there are multiple minimal values for `f(x)`, then the first one will be returned.

If `domain` is a `NamedTuple` or dictionary-like `AbstractDict` L,
the function is applied to its values.  The returned value is `L[k]` for some key `k`
such that `f(L[k])` is minimal.

##### Examples:
```julia
julia> NaNMath.argmin(abs, [1., -1., -2., 2., NaN])
1.0

julia> NaNMath.argmin(identity, [7, 1, 1, NaN])
1.0
```

julia> NaNMath.argmin(exp,Dict("x" => 1.0, "y" => -1.2, "z" => NaN))
-1.2

───────────────────────────────────────────────────────────

    NaNMath.argmin(itr) -> key

##### Args:
* `itr`: A non-empty iterable of floating point numbers.

##### Returns:
* Returns the index or key of the minimal element in `itr`, ignoring NaN's.
  If there are multiple minimal elements, then the first one will be returned.

If `itr` is a `NamedTuple` or dictionary-like `AbstractDict` L, the returned index is a key `k`,
such that `f(L[k])` is minimal.

##### Examples:
```julia
julia> NaNMath.argmin([7, 1, 1, NaN])
2

julia> NaNMath.argmin([1.0 2; 3 NaN])
CartesianIndex(1, 1)

julia> NaNMath.argmin(Dict("x" => 1.0, "y" => -1.2, "z" => NaN))
"y"
```
"""
function argmin end
argmin(f,x) = getindex(x,findmin(f,x)[2])
argmin(x) = findmin(identity,x)[2]

"""
    NaNMath.argmax(f, domain) -> x

##### Args:
* `f`: A function applied to the elements of `domain`; 
  defaults to `identity` when `domain` is the only argument.
* `domain`: A non-empty collection of floating point numbers such that
  `f` is defined on elements of `domain`.

##### Returns:
* Returns a value `x` in `domain`, ignoring NaN's, for which `f(x)` is maximized.
  If there are multiple maximal values for `f(x)`, then the first one will be returned.

If `domain` is a `NamedTuple` or dictionary-like `AbstractDict` L,
the function is applied to its values.  The returned value is `L[k]` for some key `k`
such that `f(L[k])` is maximal.

##### Examples:
```julia
julia> NaNMath.argmax(abs, [1., -1., -2., NaN])
-2.0

julia> NaNMath.argmax(identity, [7, 1, 1, NaN])
7.0
```

───────────────────────────────────────────────────────────

    NaNMath.argmax(itr) -> key

##### Args:
* `itr`: A non-empty iterable of floating point numbers.

##### Returns:
* Returns the index or key of the maximal element in `itr`, ignoring NaN's.
  If there are multiple maximal elements, then the first one will be returned.

If `itr` is a `NamedTuple` or dictionary-like `AbstractDict` L, the returned index is a key `k`,
such that `f(L[k])` is maximal.


##### Examples:
```julia
julia> NaNMath.argmax([7, 1, 1, NaN])
1

julia> NaNMath.argmax([1.0 2; 3 NaN])
CartesianIndex(2, 1)

julia> NaNMath.argmax(Dict("x" => 1.0, "y" => -1.2, "z" => NaN))
"x"
```
"""
function argmax end
argmax(x) = findmax(identity,x)[2]
argmax(f,x) = getindex(x,findmax(f,x)[2])

end
