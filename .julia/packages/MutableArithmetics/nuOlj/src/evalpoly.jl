# Copyright (c) 2023 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

# `Base.evalpoly` MA operation implemented using Horner's scheme. For
# real-valued polynomials with real coefficients.

function promote_operation(
    ::typeof(Base.evalpoly),
    ::Type{F},
    ::Type{<:Union{Tuple,AbstractVector}},
) where {F<:Real}
    return F
end

function operate_to!(
    val::F,
    ::typeof(Base.evalpoly),
    x::F,
    coefs::Union{Tuple,AbstractVector},
) where {F<:Real}
    operate!(zero, val)

    # An empty collection of coefficients is interpreted as the zero polynomial.
    isempty(coefs) && return val

    operate!(+, val, last(coefs))

    for i in eachindex(coefs)[(end-1):-1:begin]
        operate!(Base.muladd, val, x, coefs[i])
    end

    return val
end
