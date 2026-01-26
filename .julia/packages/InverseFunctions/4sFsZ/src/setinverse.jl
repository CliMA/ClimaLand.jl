# This file is a part of InverseFunctions.jl, licensed under the MIT License (MIT).


"""
    struct FunctionWithInverse{F,InvF} <: Function

A function with an inverse.

Do not construct directly, use [`setinverse(f, invf)`](@ref) instead.
"""
struct FunctionWithInverse{F,InvF} <: Function
    f::F
    invf::InvF
end
FunctionWithInverse(::Type{F}, invf::InvF) where {F,InvF} = FunctionWithInverse{Type{F},InvF}(F,invf)
FunctionWithInverse(f::F, ::Type{InvF}) where {F,InvF} = FunctionWithInverse{F,Type{InvF}}(f,InvF)
FunctionWithInverse(::Type{F}, ::Type{InvF}) where {F,InvF} = FunctionWithInverse{Type{F},Type{InvF}}(F,InvF)

(f::FunctionWithInverse)(x) = f.f(x)

inverse(f::FunctionWithInverse) = setinverse(f.invf, f.f)


"""
    setinverse(f, invf)

Return a function that behaves like `f` and uses `invf` as its inverse.

Useful in cases where no inverse is defined for `f` or to set an inverse that
is only valid within a given context, e.g. only for a limited argument
range that is guaranteed by the use case but not in general.

For example, `asin` is not a valid inverse of `sin` for arbitrary arguments
of `sin`, but can be a valid inverse if the use case guarantees that the
argument of `sin` will always be within `-π` and `π`:

```jldoctest
julia> foo = setinverse(sin, asin);

julia> x = π/3;

julia> foo(x) == sin(x)
true

julia> inverse(foo)(foo(x)) ≈ x
true

julia> inverse(foo) === setinverse(asin, sin)
true
```
"""
setinverse(f, invf) = FunctionWithInverse(_unwrap_f(f), _unwrap_f(invf))
export setinverse

_unwrap_f(f) = f
_unwrap_f(f::FunctionWithInverse) = f.f
