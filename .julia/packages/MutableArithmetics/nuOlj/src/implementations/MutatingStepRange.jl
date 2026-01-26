# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

# This file contains methods to implement the MutableArithmetics API for a new
# MutatingStepRange type. This feature was originally designed because of a
# discussion in the julialang/julia issue #39008:
# https://github.com/JuliaLang/julia/issues/39008

struct MutatingStepRange{T,S} <: OrdinalRange{T,S}
    start::T
    step::S
    stop::T
    function MutatingStepRange{T,S}(start::T, step::S, stop::T) where {T,S}
        return new(start, step, Base.steprange_last(start, step, stop))
    end
end

function MutatingStepRange(start::T, step::S, stop::T) where {T,S}
    return MutatingStepRange{T,S}(start, step, stop)
end

Base.step(r::MutatingStepRange) = r.step

function Base.unsafe_length(r::MutatingStepRange)
    n = Integer(div((r.stop - r.start) + r.step, r.step))
    return isempty(r) ? zero(n) : n
end

Base.length(r::MutatingStepRange) = Base.unsafe_length(r)

function Base.isempty(r::MutatingStepRange)
    return (r.start != r.stop) & ((r.step > zero(r.step)) != (r.stop > r.start))
end

function Base.iterate(r::MutatingStepRange)
    if isempty(r)
        return
    end
    state = copy_if_mutable(first(r))
    return (state, state)
end

function Base.iterate(r::MutatingStepRange{T}, i) where {T}
    if i == last(r)
        return
    end
    next = convert(T, operate!!(+, i, step(r)))
    return (next, next)
end
