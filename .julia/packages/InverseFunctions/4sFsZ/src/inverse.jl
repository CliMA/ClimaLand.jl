# This file is a part of InverseFunctions.jl, licensed under the MIT License (MIT).

"""
    inverse(f)

Return the inverse of function `f`.

`inverse` supports mapped and broadcasted functions (via
`Base.Broadcast.BroadcastFunction` or `Base.Fix1`) and function composition
(requires Julia >= 1.6).

# Examples

```jldoctest
julia> foo(x) = inv(exp(-x) + 1);

julia> inv_foo(y) = log(y / (1 - y));

julia> InverseFunctions.inverse(::typeof(foo)) = inv_foo;

julia> InverseFunctions.inverse(::typeof(inv_foo)) = foo;

julia> x = 4.2;

julia> inverse(foo)(foo(x)) ≈ x
true

julia> inverse(inverse(foo)) === foo
true

julia> broadcast_foo = VERSION >= v"1.6" ? Base.Broadcast.BroadcastFunction(foo) : Base.Fix1(broadcast, foo);

julia> X = rand(10);

julia> inverse(broadcast_foo)(broadcast_foo(X)) ≈ X
true

julia> bar = log ∘ foo;

julia> VERSION < v"1.6" || inverse(bar)(bar(x)) ≈ x
true
```

# Implementation

Implementations of `inverse(::typeof(f))` have to satisfy

* `inverse(f)(f(x)) ≈ x` for all `x` in the domain of `f`, and
* `inverse(inverse(f))` is defined and `inverse(inverse(f))(x) ≈ f(x)` for all `x` in the domain of `f`.

You can check your implementation with [`InverseFunctions.test_inverse`](@ref).
"""
inverse(f)
export inverse



"""
    struct NoInverse{F}

An instance `NoInverse(f)` signifies that `inverse(f)` is not defined.
"""
struct NoInverse{F}
    f::F
end
export NoInverse

NoInverse(::Type{F}) where F = NoInverse{Type{F}}(F)

(f::NoInverse)(x) = error("inverse of ", f.f, " is not defined")

inverse(f) = NoInverse(f)

inverse(f::NoInverse) = f.f



inverse(::typeof(inverse)) = inverse

@static if VERSION >= v"1.6"
    function inverse(f::Base.ComposedFunction)
        inv_inner = inverse(f.inner)
        inv_outer = inverse(f.outer)
        if inv_inner isa NoInverse || inv_outer isa NoInverse
            NoInverse(f)
        else
            Base.ComposedFunction(inv_inner, inv_outer)
        end
    end

    function inverse(bf::Base.Broadcast.BroadcastFunction)
        inv_f_kernel = inverse(bf.f)
        if inv_f_kernel isa NoInverse
            NoInverse(bf)
        else
            Base.Broadcast.BroadcastFunction(inv_f_kernel)
        end
    end
end

function inverse(mapped_f::Base.Fix1{<:Union{typeof(map),typeof(broadcast)}})
    inv_f_kernel = inverse(mapped_f.x)
    if inv_f_kernel isa NoInverse
        NoInverse(mapped_f)
    else
        Base.Fix1(mapped_f.f, inv_f_kernel)
    end
end

inverse(::typeof(identity)) = identity
inverse(::typeof(inv)) = inv
inverse(::typeof(adjoint)) = adjoint
inverse(::typeof(transpose)) = transpose
inverse(::typeof(conj)) = conj

inverse(::typeof(!)) = !
inverse(::typeof(+)) = +
inverse(::typeof(-)) = -

inverse(f::Base.Fix1{typeof(+)}) = Base.Fix2(-, f.x)
inverse(f::Base.Fix2{typeof(+)}) = Base.Fix2(-, f.x)
inverse(f::Base.Fix1{typeof(-)}) = Base.Fix1(-, f.x)
inverse(f::Base.Fix2{typeof(-)}) = Base.Fix1(+, f.x)
inverse(f::Base.Fix1{typeof(*)}) = iszero(f.x) ? throw(DomainError(f.x, "Cannot invert multiplication by zero")) : Base.Fix1(\, f.x)
inverse(f::Base.Fix2{typeof(*)}) = iszero(f.x) ? throw(DomainError(f.x, "Cannot invert multiplication by zero")) : Base.Fix2(/, f.x)
inverse(f::Base.Fix1{typeof(/)}) = Base.Fix2(\, f.x)
inverse(f::Base.Fix2{typeof(/)}) = Base.Fix2(*, f.x)
inverse(f::Base.Fix1{typeof(\)}) = Base.Fix1(*, f.x)
inverse(f::Base.Fix2{typeof(\)}) = Base.Fix1(/, f.x)

inverse(::typeof(deg2rad)) = rad2deg
inverse(::typeof(rad2deg)) = deg2rad

inverse(::typeof(exp)) = log
inverse(::typeof(log)) = exp

inverse(::typeof(exp2)) = log2
inverse(::typeof(log2)) = exp2

inverse(::typeof(exp10)) = log10
inverse(::typeof(log10)) = exp10

inverse(::typeof(expm1)) = log1p
inverse(::typeof(log1p)) = expm1


inverse(::typeof(sinh)) = asinh
inverse(::typeof(tanh)) = atanh
inverse(::typeof(coth)) = acoth
inverse(::typeof(csch)) = acsch

inverse(::typeof(asinh)) = sinh
inverse(::typeof(atanh)) = tanh
inverse(::typeof(acoth)) = coth
inverse(::typeof(acsch)) = csch


inverse(::typeof(sqrt)) = square
inverse(::typeof(square)) = sqrt

inverse(::typeof(cbrt)) = Base.Fix2(^, 3)

inverse(f::Base.Fix2{typeof(^)}) = iszero(f.x) ? throw(DomainError(f.x, "Cannot invert x^$(f.x)")) : Base.Fix2(invpow_arg2, f.x)
inverse(f::Base.Fix2{typeof(^), <:Integer}) = isodd(f.x) ? Base.Fix2(invpow_arg2, f.x) : throw(DomainError(f.x, "Cannot invert x^$(f.x)"))
inverse(f::Base.Fix2{typeof(invpow_arg2)}) = Base.Fix2(^, f.x)

inverse(f::Base.Fix1{typeof(^), <:Real}) = f.x > zero(f.x) ? Base.Fix1(invpow_arg1, f.x) : throw(DomainError(f.x, "Cannot invert $(f.x)^x"))
inverse(f::Base.Fix1{typeof(^)}) = Base.Fix1(invpow_arg1, f.x)
inverse(f::Base.Fix1{typeof(invpow_arg1)}) = Base.Fix1(^, f.x)
inverse(f::Base.Fix1{typeof(log)}) = isone(f.x) ? throw(DomainError(f.x, "Cannot invert log($(f.x), x)")) : Base.Fix1(invlog_arg1, f.x)
inverse(f::Base.Fix1{typeof(invlog_arg1)}) = Base.Fix1(log, f.x)

inverse(f::Base.Fix2{typeof(log)}) = isone(f.x) ? throw(DomainError(f.x, "Cannot invert log(x, $(f.x))")) : Base.Fix2(invlog_arg2, f.x)
inverse(f::Base.Fix2{typeof(invlog_arg2)}) = Base.Fix2(log, f.x)


inverse(f::Base.Fix2{typeof(divrem)}) = Base.Fix2(invdivrem, f.x)
inverse(f::Base.Fix2{typeof(invdivrem)}) = Base.Fix2(divrem, f.x)

inverse(f::Base.Fix2{typeof(fldmod)}) = Base.Fix2(invfldmod, f.x)
inverse(f::Base.Fix2{typeof(invfldmod)}) = Base.Fix2(fldmod, f.x)

inverse(::typeof(reim)) = Base.splat(complex)
inverse(::typeof(Base.splat(complex))) = reim

inverse(::typeof(reverse)) = reverse
