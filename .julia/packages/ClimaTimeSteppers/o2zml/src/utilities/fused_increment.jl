"""
    fused_increment(u, dt, sc::SparseCoeffs, tend, ::Val{i}) where {i}

Returns a broadcasted object in the form

    `u + ∑ⱼ dt scⱼ  tendⱼ for j in 1:i`
    or
    `u + ∑ⱼ dt scᵢⱼ tendⱼ for j in 1:(i-1)`

depending on the dimensions of the coefficients in `sc`. The
broadcasted object drops zero coefficients from the expression
at compile-time using the mask `mask` (made from `SparseCoeffs(coeffs)`).

Returns `u` when `i ≤ 1` or if the mask is all true (all coefficients are zero).


# Example loops that this is fusing:

For 2D coefficients case
```julia
for j in 1:(i - 1)
    iszero(a_imp[i, j]) && continue
    @. U += dt * a_imp[i, j] * T_imp[j]
end

For 1D coefficients case
```julia
for j in 1:s
    iszero(b_exp[j]) && continue
    @. temp += dt * b_exp[j] * T_lim[j]
end
```
"""
function fused_increment end

@inline fused_increment(u, dt, sc, tend, v) = fused_increment(u, dt, sc, tend, v, get_S(sc))

# =================================================== ij case (S::NTuple{2})
# recursion: _rfused_increment_j always returns a Tuple
@inline _rfused_increment_ij(js::Tuple{}, i, u, dt, sc::SparseCoeffs, tend) = ()

@inline _rfused_increment_ij(js::Tuple{Int}, i, u, dt, sc::T, tend) where {T <: SparseCoeffs} =
    _rfused_increment_ij(js[1], i, u, dt, sc, tend)

@inline _rfused_increment_ij(j::Int, i, u, dt, sc::T, tend) where {T <: SparseCoeffs} =
    if zero_coeff(T, i, j)
        ()
    else
        (Base.Broadcast.broadcasted(*, dt * sc[i, j], tend[j]),)
    end

@inline _rfused_increment_ij(js::Tuple, i, u, dt, sc::T, tend) where {T <: SparseCoeffs} =
    (_rfused_increment_ij(first(js), i, u, dt, sc, tend)..., _rfused_increment_ij(Base.tail(js), i, u, dt, sc, tend)...)

# top-level function
@inline function _fused_increment_ij(js::Tuple, i, u, dt, sc::T, tend) where {T <: SparseCoeffs}
    return if all(j -> zero_coeff(T, i, j), js)
        u
    else
        Base.Broadcast.broadcasted(
            +,
            u,
            _rfused_increment_ij(js, i, u, dt, sc, tend)..., # recurse...
        )
    end
end

# top-level function, if the tuple is empty, just return u
@inline _fused_increment_ij(js::Tuple{}, i, u, dt, sc::SparseCoeffs, tend) = u

# wrapper ij case (S::NTuple{2})
@inline fused_increment(u, dt, sc::SparseCoeffs, tend, ::Val{i}, ::NTuple{2}) where {i} =
    _fused_increment_ij(ntuple(j -> j, Val(i - 1)), i, u, dt, sc, tend)

# =================================================== j case (S::NTuple{1})
# recursion: _rfused_increment_j always returns a Tuple
@inline _rfused_increment_j(js::Tuple{}, u, dt, sc::SparseCoeffs, tend) = ()

@inline _rfused_increment_j(js::Tuple{Int}, u, dt, sc::T, tend) where {T <: SparseCoeffs} =
    _rfused_increment_j(js[1], u, dt, sc, tend)

@inline _rfused_increment_j(j::Int, u, dt, sc::T, tend) where {T <: SparseCoeffs} =
    if zero_coeff(T, j)
        ()
    else
        (Base.Broadcast.broadcasted(*, dt * sc[j], tend[j]),)
    end

@inline _rfused_increment_j(js::Tuple, u, dt, sc::T, tend) where {T <: SparseCoeffs} =
    (_rfused_increment_j(first(js), u, dt, sc, tend)..., _rfused_increment_j(Base.tail(js), u, dt, sc, tend)...)

# top-level function
@inline function _fused_increment_j(js::Tuple, u, dt, sc::T, tend) where {T <: SparseCoeffs}
    return if all(j -> zero_coeff(T, j), js)
        u
    else
        Base.Broadcast.broadcasted(
            +,
            u,
            _rfused_increment_j(js, u, dt, sc, tend)..., # recurse...
        )
    end
end

# top-level function, if the tuple is empty, just return u
@inline _fused_increment_j(js::Tuple{}, u, dt, sc::SparseCoeffs, tend) = u

# wrapper j case (S::NTuple{1})
@inline fused_increment(u, dt, sc::SparseCoeffs, tend, ::Val{s}, ::NTuple{1}) where {s} =
    _fused_increment_j(ntuple(i -> i, s), u, dt, sc, tend)

"""
    fused_increment!(u, dt, sc, tend, v)

Calls [`fused_increment`](@ref) and materializes
a broadcast expression in the form:

    `@. u += ∑ⱼ dt scⱼ tendⱼ`

In the edge case (coeffs are zero, `j` range is empty),
this lowers to `nothing` (no-op)
"""
@inline function fused_increment!(u, dt, sc, tend, v)
    bc = fused_increment(u, float(dt), sc, tend, v)
    if bc isa Base.Broadcast.Broadcasted # Only material if not trivial assignment
        Base.Broadcast.materialize!(u, bc)
    end
    nothing
end

"""
    assign_fused_increment!(U, u, dt, sc, tend, v)

Calls [`fused_increment`](@ref) and materializes
a broadcast expression in the form:

    `@. u += ∑ⱼ dt scⱼ tendⱼ`

In the edge case (coeffs are zero, `j` range is empty),
this lowers to

    `@. U = u`
"""
@inline function assign_fused_increment!(U, u, dt, sc, tend, v)
    bc = fused_increment(u, float(dt), sc, tend, v)
    Base.Broadcast.materialize!(U, bc)
    return nothing
end
