# Because there exists the `FiniteDiff.jl` package with quite different purposes,
# this module is not expected to be reexported.
module FiniteDiff

using ImageCore
using ImageCore.MappedArrays: of_eltype

"""
Although stored as an array, image can also be viewed as a function from discrete grid space
Zᴺ to continuous space R if it is gray image, to C if it is complex-valued image
(MRI rawdata), to Rᴺ if it is colorant image, etc.
This module provides the discrete
version of gradient-related operators by viewing image arrays as functions.

This module provides:

- forward/backward difference [`fdiff`](@ref) are the Images-flavor of `Base.diff`
- gradient operator [`fgradient`](@ref) and its adjoint via keyword `adjoint=true`.
- divergence operator [`fdiv`](@ref) computes the sum of discrete derivatives of vector
  fields.
- laplacian operator [`flaplacian`](@ref) is the divergence of the gradient fields.

Every function in this module has its in-place version.
"""
FiniteDiff

export
    fdiff, fdiff!,
    fdiv, fdiv!,
    fgradient, fgradient!,
    flaplacian, flaplacian!


"""
    fdiff(A::AbstractArray; dims::Int, rev=false, boundary=:periodic)

A one-dimension finite difference operator on array `A`. Unlike `Base.diff`, this function doesn't
shrink the array size.

Take vector as an example, it computes `(A[2]-A[1], A[3]-A[2], ..., A[1]-A[end])`.

# Keywords

- `rev::Bool`
  If `rev==true`, then it computes the backward difference
  `(A[end]-A[1], A[1]-A[2], ..., A[end-1]-A[end])`.
- `boundary`
  By default it computes periodically in the boundary, i.e., `:periodic`.
  In some cases, one can fill zero values with `boundary=:zero`.

# Examples

```jldoctest; setup=:(using ImageBase.FiniteDiff: fdiff)
julia> A = [2 4 8; 3 9 27; 4 16 64]
3×3 $(Matrix{Int}):
 2   4   8
 3   9  27
 4  16  64

julia> diff(A, dims=2) # this function exists in Base
3×2 $(Matrix{Int}):
  2   4
  6  18
 12  48

julia> fdiff(A, dims=2)
3×3 $(Matrix{Int}):
  2   4   -6
  6  18  -24
 12  48  -60

julia> fdiff(A, dims=2, rev=true) # reverse diff
3×3 $(Matrix{Int}):
  -6   2   4
 -24   6  18
 -60  12  48

julia> fdiff(A, dims=2, boundary=:zero) # fill boundary with zeros
3×3 $(Matrix{Int}):
  2   4  0
  6  18  0
 12  48  0
```

See also [`fdiff!`](@ref) for the in-place version.
"""
fdiff(A::AbstractArray; kwargs...) = fdiff!(similar(A, maybe_floattype(eltype(A))), A; kwargs...)

"""
    fdiff!(dst::AbstractArray, src::AbstractArray; dims::Int, rev=false, boundary=:periodic)

The in-place version of [`fdiff`](@ref)
"""
function fdiff!(dst::AbstractArray, src::AbstractArray;
        dims=_fdiff_default_dims(src),
        rev=false,
        boundary::Symbol=:periodic)
    isnothing(dims) && throw(UndefKeywordError(:dims))
    axes(dst) == axes(src) || throw(ArgumentError("axes of all input arrays should be equal. Instead they are $(axes(dst)) and $(axes(src))."))
    N = ndims(src)
    1 <= dims <= N || throw(ArgumentError("dimension $dims out of range (1:$N)"))

    src = of_eltype(maybe_floattype(eltype(dst)), src)
    r = axes(src)
    r0 = ntuple(i -> i == dims ? UnitRange(first(r[i]), last(r[i]) - 1) : UnitRange(r[i]), N)
    r1 = ntuple(i -> i == dims ? UnitRange(first(r[i])+1, last(r[i])) : UnitRange(r[i]), N)

    d0 = ntuple(i -> i == dims ? UnitRange(last(r[i]), last(r[i])) : UnitRange(r[i]), N)
    d1 = ntuple(i -> i == dims ? UnitRange(first(r[i]), first(r[i])) : UnitRange(r[i]), N)

    if rev
        dst[r1...] .= view(src, r1...) .- view(src, r0...)
        if boundary == :periodic
            dst[d1...] .= view(src, d1...) .- view(src, d0...)
        elseif boundary == :zero
            dst[d1...] .= zero(eltype(dst))
        else
            throw(ArgumentError("Wrong boundary condition $boundary"))
        end
    else
        dst[r0...] .= view(src, r1...) .- view(src, r0...)
        if boundary == :periodic
            dst[d0...] .= view(src, d1...) .- view(src, d0...)
        elseif boundary == :zero
            dst[d0...] .= zero(eltype(dst))
        else
            throw(ArgumentError("Wrong boundary condition $boundary"))
        end
    end

    return dst
end

_fdiff_default_dims(A) = nothing
_fdiff_default_dims(A::AbstractVector) = 1

maybe_floattype(::Type{T}) where T = T
maybe_floattype(::Type{T}) where T<:FixedPoint = floattype(T)
maybe_floattype(::Type{CT}) where CT<:Color = base_color_type(CT){maybe_floattype(eltype(CT))}


"""
    fdiv(Vs...)

Compute the divergence of vector fields `Vs`.

See also [`fdiv!`](@ref) for the in-place version.
"""
function fdiv(V₁::AbstractArray{T}, Vs::AbstractArray{T}...) where T
    fdiv!(similar(V₁, maybe_floattype(T)), V₁, Vs...)
end
fdiv(Vs::Tuple) = fdiv(Vs...)

"""
    fdiv!(out, Vs...)

The in-place version of divergence operator [`fdiv`](@ref).
"""
function fdiv!(out::AbstractArray, V₁::AbstractArray, Vs::AbstractArray...)
    # This is the optimized version of `sum(v->fgradient(v; ajoint=true), (V₁, Vs...))`
    # by removing unnecessary memory allocations.
    all(v->axes(v) == axes(out), (V₁, Vs...)) || throw(ArgumentError("All axes of vector fields Vs and X should be the same."))

    # TODO(johnnychen94): for better performance, we can eliminate this `tmp` allocation by fusing multiple `fdiff` in the inner loop.
    out .= fdiff(V₁; dims=1, rev=true, boundary=:periodic)
    tmp = similar(out)
    for i in 1:length(Vs)
        out .+= fdiff!(tmp, Vs[i]; dims=i+1, rev=true, boundary=:periodic)
    end
    return out
end
fdiv!(out::AbstractArray, Vs::Tuple) = fdiv!(out, Vs...)

"""
    flaplacian(X::AbstractArray)

The discrete laplacian operator, i.e., the divergence of the gradient fields of `X`.

See also [`flaplacian!`](@ref) for the in-place version.
"""
flaplacian(X::AbstractArray) = flaplacian!(similar(X, maybe_floattype(eltype(X))), X)

"""
    flaplacian!(out, X)
    flaplacian!(out, ∇X::Tuple, X)

The in-place version of the laplacian operator [`flaplacian`](@ref).

!!! tip Avoiding allocations
    The two-argument method will allocate memory to store the intermediate
    gradient fields `∇X`. If you call this repeatedly with images of consistent size and type,
    consider using the three-argument form with pre-allocated memory for `∇X`,
    which will eliminate allocation by this function.
"""
flaplacian!(out, X::AbstractArray) = fdiv!(out, fgradient(X))
flaplacian!(out, ∇X::Tuple, X::AbstractArray) = fdiv!(out, fgradient!(∇X, X))


"""
    fgradient(X::AbstractArray; adjoint=false, boundary=:periodic) -> (∂₁X, ∂₂X, ..., ∂ₙX)

Computes the gradient fields of `X`. If `adjoint==true` then it computes the adjoint gradient
fields.

Each gradient vector is computed as forward difference along specific dimension, e.g.,
[`∂ᵢX = fdiff(X, dims=i)`](@ref fdiff). The `boundary` keyword is passed to `fdiff` to
specify the behavior on boundary.

Mathematically, the adjoint operator ∂ᵢ' of ∂ᵢ is defined as `<∂ᵢu, v> := <u, ∂ᵢ'v>`, where
`<a, b>` is the inner product of vector `a` and `b`.

See also the in-place version [`fgradient!(X)`](@ref) to reuse the allocated memory.
"""
function fgradient(X::AbstractArray{T,N}; adjoint::Bool=false, boundary=:periodic) where {T,N}
    fgradient!(ntuple(i->similar(X, maybe_floattype(T)), N), X; adjoint=adjoint, boundary=boundary)
end

"""
    fgradient!(∇X::Tuple, X::AbstractArray; adjoint=false, boundary=:periodic)

The in-place version of (adjoint) gradient operator [`fgradient`](@ref).

The input `∇X = (∂₁X, ∂₂X, ..., ∂ₙX)` is a tuple of arrays that are similar to `X`, i.e.,
`eltype(∂ᵢX) == eltype(X)` and `axes(∂ᵢX)  == axes(X)` for all `i`.
"""
function fgradient!(∇X::NTuple{N, <:AbstractArray}, X; adjoint::Bool=false, boundary=:periodic) where N
    all(v->axes(v) == axes(X), ∇X) || throw(ArgumentError("All axes of vector fields ∇X and X should be the same."))
    for i in 1:N
        if adjoint
            # the negative adjoint of gradient operator for forward difference is the backward difference
            # see also
            # Getreuer, Pascal. "Rudin-Osher-Fatemi total variation denoising using split Bregman." _Image Processing On Line_ 2 (2012): 74-95.
            fdiff!(∇X[i], X, dims=i, rev=true, boundary=boundary)
            # TODO(johnnychen94): ideally we can get avoid flipping the signs for better performance.
            @. ∇X[i] = -∇X[i]
        else
            fdiff!(∇X[i], X, dims=i, boundary=boundary)
        end
    end
    return ∇X
end



if VERSION < v"1.1"
    isnothing(x) = x === nothing
end

end # module
