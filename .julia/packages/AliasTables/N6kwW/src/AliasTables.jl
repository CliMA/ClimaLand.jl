module AliasTables

using Random, PtrArrays

export AliasTable
VERSION >= v"1.11.0-DEV.469" && eval(Meta.parse("public sample, probabilities, set_weights!"))

const Memory = isdefined(Base, :Memory) ? Base.Memory : Vector # VERSION <= 1.10

if isdefined(Base, :top_set_bit)
    const top_set_bit = Base.top_set_bit
else
    top_set_bit(x::Integer) = top_set_bit(UInt64(x)) # VERSION <= 1.9
    top_set_bit(x::Base.BitInteger) = 8sizeof(x) - leading_zeros(x)
end

if isdefined(Base, :require_one_based_indexing) # VERSION == 1.0
    const require_one_based_indexing = Base.require_one_based_indexing
else
    require_one_based_indexing(A...) = !Base.has_offset_axes(A...) || throw(ArgumentError("offset arrays are not supported but got an array with index other than 1"))
end

"""
    AliasTable{T<:Unsigned=UInt64, I<:Integer=Int}(weights::AbstractVector{<:Real})

An efficient data structure for sampling from a discrete distribution.

Maps every value representable by `T` to a value of type `I` in `eachindex(wights)` such
that the number of values maped to a given index of `weights` is proportional to the value
at that index.

The mapping can be accessed directly via
[`AliasTables.sample(x::T, at::AliasTable{T, I})`](@ref AliasTables.sample)
or indirectly via the `Random` API: `rand(at)`, `rand(rng, at)`, `rand(at, dims...)`, etc.

See also [`AliasTables.set_weights!`](@ref)

# Example

```jldoctest; filter=[r" [1-3]"]
julia> at = AliasTable([1, 3, 1])
AliasTable([0x3333333333333334, 0x9999999999999999, 0x3333333333333333])

julia> rand(at, 5)
5-element Vector{Int64}:
 2
 3
 2
 2
 3
```
"""
struct AliasTable{T <: Unsigned, I <: Integer}
    mask::T
    probability_alias::Memory{Tuple{T, I}}
    length::I

    function AliasTable{T, I}(weights::AbstractVector{<:Real}; _normalize=true) where {T <: Unsigned, I <: Integer}
        length(weights) > typemax(I) && throw(ArgumentError(
            "length(weights) must be less than typemax(I). Got $(length(weights)) and $(typemax(I)), respectively."))
        shift = top_set_bit(length(weights) - 1)
        probability_alias = Memory{Tuple{T, I}}(undef, 1 << shift)
        mask = (one(T) << max(8sizeof(T) - shift, 0)) - one(T)
        at = new{T, I}(mask, probability_alias, length(weights))
        set_weights!(at, weights, _normalize=_normalize)
    end
end

AliasTable(weights::AbstractVector{<:Real}; _normalize=true) = AliasTable{UInt64, Int}(weights; _normalize=_normalize)
AliasTable{T}(weights::AbstractVector{<:Real}; _normalize=true) where T <: Unsigned = AliasTable{T, Int}(weights; _normalize=_normalize)

# The constructors (i.e. set_weights! and all its special cases) are responsible for
# ensuring that `at1.probability_alias == at2.probability_alias` whenever `at1` and
# `at2` have the same weights according to `AliasTables.probabilities(float, at)`

"""
    set_weights!(at::AliasTable, weights::AbstractVector{<:Real})

Set the weights of `at` to `weights` and return `at`.

Does not perform GC managed allocations.

# Example
```jldoctest
julia> at = AliasTable([1, 3, 1])
AliasTable([0x3333333333333334, 0x9999999999999999, 0x3333333333333333])

julia> at === AliasTables.set_weights!(at, [1, 2, 1])
true

julia> at
AliasTable([0x4000000000000000, 0x8000000000000000, 0x4000000000000000])
```
"""
function set_weights!(at::AliasTable{T}, weights::AbstractVector{<:Real}; _normalize=true) where T
    require_one_based_indexing(weights)
    length(weights) == length(at) || throw(DimensionMismatch("length(weights) must equal length(at). Got $(length(weights)) and $(length(at)), respectively."))
    probability_alias = at.probability_alias
    if _normalize
        (is_constant, sm) = checked_sum(weights)
        if is_constant
            _constant_alias_table!(probability_alias, sm, length(weights))
        elseif sm-true == typemax(T) # pre-normalized
            _alias_table!(probability_alias, weights)
        elseif iszero(sm)
            throw(ArgumentError("sum(weights) overflows, but just barely"))
        else
            norm = malloc(T, length(weights))
            try
                normalize_to_uint!(norm, weights, sm)
                _alias_table!(probability_alias, norm)
            finally
                free(norm)
            end
        end
    else
        throw_on_negatives(weights)
        _alias_table!(probability_alias, weights)
    end
    at
end

function _constant_alias_table!(probability_alias::Memory{Tuple{T, I}}, index, length) where {T, I}
    amount_redirected = one(T) << (8*sizeof(T) - 1) # typemax(T)+1 / 2, typically more than points per cell!
    for i in eachindex(probability_alias)
        probability_alias[i] = (amount_redirected, (index-i)%I)
    end
end

function _lookup_alias_table!(probability_alias::Memory{Tuple{T, I}}, weights, mtz) where {T, I}
    amount_redirected = one(T) << (8*sizeof(T) - 1) # typemax(T)+1 / 2, typically more than points per cell!
    len = 1 << (8sizeof(T) - mtz)
    len == 0 && throw(ArgumentError("Lookup table longer than typemax(Int)"))
    len > length(probability_alias) && throw(ArgumentError("Lookup table longer than length(probablity_alias)"))
    j = 1
    for (i, w) in enumerate(weights)
        j2 = j + Int(w >> mtz)
        j2-1 > len && throw(ArgumentError("sum(weights) is too high"))
        for k in j:j2-1
            probability_alias[k] = (amount_redirected, (i-k)%I)
        end
        j = j2
    end
    j <= len && throw(ArgumentError("sum(weights) is too low"))
    for j in len+1:length(probability_alias)
        (p,a) = probability_alias[j-len]
        probability_alias[j] = p, (a-len)%I
    end
end

function minimum_trailing_zeros(x) # minimum(trailing_zeros, x), but faster.
    v, i = iterate(x)
    orv = v
    vi = v, i
    while iseven(orv)
        vi = iterate(x, i)
        vi === nothing && break
        v, i = vi
        orv |= v
    end
    trailing_zeros(orv)
end

function throw_on_negatives(weights)
    for w in weights
        w < 0 && throw(ArgumentError("found negative weight $w"))
    end
end
function get_only_nonzero(weights)
    only_nonzero = -1
    for (i, w) in enumerate(weights)
        if w > 0
            if only_nonzero == -1
                only_nonzero = i
            else
                only_nonzero = -2
                break
            end
        end
    end
    only_nonzero == -1 && throw(ArgumentError("all weights are zero"))
    only_nonzero
end

function _alias_table!(probability_alias::Memory{Tuple{T, I}}, weights) where {T, I}
    onz = get_only_nonzero(weights)
    onz == -2 || return _constant_alias_table!(probability_alias, onz, length(weights))

    bitshift = top_set_bit(length(weights) - 1)
    len = 1 << bitshift # next_or_eq_power_of_two(length(weights))
    @assert len == length(probability_alias)
    points_per_cell = one(T) << (8*sizeof(T) - bitshift) # typemax(T)+1 / len

    # The reason this optimiation exists is to prevent when points_per_cell == 0.
    # The reason it is applied so aggressively is so aggressively is so that
    # AliasTables of different bitwidths are structurally similar for the same
    # weights. This is important for equality and hashing.
    mtz = minimum_trailing_zeros(weights)
    lookup_table_bits = 8sizeof(T) - mtz
    lookup_table_bits < bitshift && return _lookup_alias_table!(probability_alias, weights, mtz)

    # @show sum(weights)
    # @show weights

    enum_weights = enumerate(weights)
    (surplus_i, surplus_desired), surplus_state = (thirsty_i, thirsty_desired), thirsty_state = iterate(enumerate(weights))
    current_i = surplus_i
    current_desired = surplus_desired

    while true
        # @show current_i, current_desired, points_per_cell
        if current_desired < points_per_cell # Surplus (strict)
            while true                                                            # Find the next thirsty cell
                ix = iterate(enum_weights, thirsty_state)
                ix === nothing && throw(ArgumentError("sum(weights) is too low")) # If there is no thirsty cell, there are more points than requsted by weights.
                (thirsty_i, thirsty_desired), thirsty_state = ix
                thirsty_desired >= points_per_cell && break
            end
            excess = points_per_cell - current_desired                            # Assign this many extra points
            probability_alias[current_i] = (excess, (thirsty_i-current_i)%I)      # To the targeted cell
            current_i = thirsty_i                                                 # Now we have to make sure that thristy cell gets exactly what it wants and no more
            current_desired = thirsty_desired - excess                            # It wants what it wants and hasn't already been transferred
        else                                 # Thirsty (loose)
            while true                                                            # Find the next surplus cell
                ix = iterate(enum_weights, surplus_state)
                ix === nothing && @goto break_outer                               # If there is no surplus cell, handle below
                (surplus_i, surplus_desired), surplus_state = ix
                surplus_desired < points_per_cell && break
            end
            excess = points_per_cell - surplus_desired                            # Assign this many extra points
            probability_alias[surplus_i] = (excess, (current_i-surplus_i)%I)      # From the cell with surplus to this cell
            current_desired -= excess                                             # We now don't want as many points (and may even no longer be thristy)
        end
    end
    @label break_outer
    # println()

    # There are no real surplus cells, but there may be synthetic surplus cells to round out
    # to a power of two. Those synthetic cells have structural weight of 0 so surplus value
    # of points_per_cell each. There may be remaining thristy cells, and the current cell is
    # thirsty.

    surplus_state_2 = length(weights)


    while true
        # @show current_i, current_desired, points_per_cell
        if current_desired < points_per_cell # Surplus (strict)
            while true                                                                # Find the next thirsty cell
                ix = iterate(enum_weights, thirsty_state)
                ix === nothing && throw(ArgumentError("sum(weights) is too low"))     # If there is no thirsty cell, there are more points than requsted by weights.
                (thirsty_i, thirsty_desired), thirsty_state = ix
                thirsty_desired >= points_per_cell && break
            end
            excess = points_per_cell - current_desired                                # Assign this many extra points
            probability_alias[current_i] = (excess, (thirsty_i-current_i)%I)          # To the targeted cell
            current_i = thirsty_i                                                     # Now we have to make sure that thristy cell gets exactly what it wants and no more
            current_desired = thirsty_desired - excess                                # It wants what it wants and hasn't already been transferred
        else                                 # Thirsty (loose)
            surplus_state_2 += true # Find the next surplus cell
            surplus_state_2 > len && break # If there is no surplus cell, handle below
            surplus_i = surplus_state_2
            excess = points_per_cell                                                  # Assign all the points
            probability_alias[surplus_i] = (points_per_cell, (current_i-surplus_i)%I) # From the synthetic cell with surplus to this cell
            current_desired -= excess                                                 # We now don't want as many points (and may even no longer be thristy)
        end
    end

    # @show probability_alias, points_per_cell, current_desired

    if points_per_cell < current_desired  # Strictly thirsty, and no surplus cells, so exceed the desired weight.
        throw(ArgumentError("sum(weights) is too high"))
    end
    probability_alias[current_i] = (0, 0)
    # Just right. There are no surplus cells left and no current surplus or thirst. All
    # that's left are future loosely thirsty cells, all of which should be a thirst of
    # exactly 0.
    while true         # Loop over all thirsty celss
        ix = iterate(enum_weights, thirsty_state)
        ix === nothing && break # Out of thirsty cells, yay!
        (thirsty_i, thirsty_desired), thirsty_state = ix
        points_per_cell < thirsty_desired && throw(ArgumentError("sum(weights) is too high")) # Strictly thirsty, with no surplus to draw from.
        points_per_cell == thirsty_desired && (probability_alias[thirsty_i] = (0, 0)) # loosely thirsty, but satisfied. Zero out the undef.
    end
end

"""
    sample(x::T, at::AliasTable{T, I}) -> I

Sample from `at` using the seed `x`.

If `x` is chosen uniformly at random from the set of all values representable by `T` then
the output will be a random sample from the distribution represented by `at`. The mapping is
deterministic and not pseudo-random so for patterned input `x` the output will be patterned
as well.

See also [`AliasTable`](@ref), [`AliasTables.probabilities`](@ref)
"""
function sample(x::T, at::AliasTable{T, I}) where {T, I}
    shift = max(top_set_bit(typemax(T)) - top_set_bit(length(at.probability_alias)) + 1, 0)
    cell = (x >> shift) + 1
    # @assert (one(T) << shift) - one(T) == at.mask
    val = x & at.mask
    prob, alias = @inbounds at.probability_alias[cell%Int] # see proof below
    (((val < prob) * alias + cell)%I)::I
end

function _sample!(xs::AbstractArray{I}, at::AliasTable{T, I}) where {T<:Base.BitUnsigned, I<:Base.BitInteger}
    @assert sizeof(I) <= sizeof(T)
    shift = max(top_set_bit(typemax(T)) - top_set_bit(length(at.probability_alias)) + 1, 0)
    # @assert (one(T) << shift) - one(T) == at.mask
    @simd ivdep for i in eachindex(xs)
        x = unsigned(xs[i])
        cell = (x >> shift) + 1
        val = x & at.mask
        prob, alias = @inbounds at.probability_alias[cell%Int] # see proof below
        xs[i] = (((val < prob) * alias + cell)%I)::I
    end
    xs
end

function Random.rand!(rng::Random.AbstractRNG,
                      dst::Array{I},
                      st::Random.SamplerTrivial{AliasTable{T, I}, I}) where {T<:Base.BitUnsigned, I<:Base.BitInteger}
    if sizeof(I) == sizeof(T)
        rand!(rng, dst)
        _sample!(dst, st.self)
    else
        for i in eachindex(dst)
            dst[i] = rand(rng, st)
        end
    end
    dst
end

#=
Loose justification that @inbounds is safe:

The worry is that `x` could be too large, so let's assume `x == typemax(T)`. When we perform
the bitshift `x >> shift`, we can decompose that into shifting by each of summands of shift,
one after another, but only dropping out of bounds bits at the end. Lookinf first at
`x >> top_set_bit(typemax(T))`, this is, in the worst case, just below 1. Now, let's assume
that `length(at.probability_alias)` is a power of 2. In this case,
`length(at.probability_alias) == 1 << (top_set_bit(length(at.probability_alias))-1)`
and we know that `cell = ([something less than 1] >> (-top_set_bit(length(at.probability_alias))+1)) + 1`,
so `cell < length(at.probability_alias) + 1`, and therefore `cell` is in bounds. Moding by
`Int` isn't going make `cell` any bigger. In the event that `length(at.probability_alias)`
is not a power of 2, the bounds we established above still hold based on the top set bit of
that length, but some lower bits are also 1s which gives us some extra room to spare.

Mathematical proof that the @inbounds is safe:

Bitshifts with unsigned left hand sides have the property that `x >> n` is always less than
or equal to the mathematical value ``x2^{-n}``, with equality if none of the bits are
shifted off the edge of the word. This holds regardless of the sign on `n`.

Consequently, ``cell-1 = x >> shift ≤ x2^{-shift} ≤ typemax(T)2^{-shift}``.

Let us abbreviate `top_set_bit` as `tsb`, and note that we always
have ``2^{tsb(x)-1} ≤ x < 2^{tsb(x)}``. This lets us write

``typemax(T)2^{-shift} < 2^{tsb(typemax(T))}2^{-shift} = 2^{tsb(typemax(T))-shift}``

Expanding `shift` we find

``tsb(typemax(T))-shift = tsb(typemax(T))-(tsb(typemax(T)) - tsb(length(at.probability_alias)) + 1)
                        = tsb(length(at.probability_alias))-1``

So we have ``2^{tsb(typemax(T))-shift} = 2^{tsb(length(at.probability_alias))-1}``.

From the innequality above, `2^{tsb(length(at.probability_alias))-1} ≤ length(at.probability_alias))`

Stringing all these together we have ``cell-1 < length(at.probability_alias)``, and
``cell-1`` is clearly non-negative, so we have ``1 ≤ cell ≤ length(at.probability_alias)``.

Assuming that `length(at.probability_alias)` is an `Int`, `cell%Int == cell`, and indexing
into `at.probability_alias[cell%Int]` is in bounds.
=#

### Random API
Random.rand(rng::Random.AbstractRNG, at::Random.SamplerTrivial{<:AliasTable{T}}) where T = sample(rand(rng, T), at.self)
Random.gentype(::Type{AliasTable{T, I}}) where {T, I} = I

### Reconstruct probabilities
"""
    probabilities(at::AliasTable{T}) -> Vector{T}

Recover the exact sampling weights from a given `AliasTable`. The returned values will
sum to one more than `typemax(T)`, unless `at` is a constant distribution (e.g.
`AliasTable([0,1,0])`), in which case the weights will sum to `typemax(T)`.

See also [`AliasTable`](@ref), [`AliasTables.sample`](@ref)

# Examples

```jldoctest
julia> at = AliasTable([1, 3, 1])
AliasTable([0x3333333333333334, 0x9999999999999999, 0x3333333333333333])

julia> AliasTables.probabilities(at)
3-element Vector{UInt64}:
 0x3333333333333334
 0x9999999999999999
 0x3333333333333333

julia> AliasTables.probabilities(AliasTable([0, 1, 0]))
3-element Vector{UInt64}:
 0x0000000000000000
 0xffffffffffffffff
 0x0000000000000000
```
"""
function probabilities(at::AliasTable{T, I}) where {T, I}
    bitshift = top_set_bit(length(at.probability_alias) - 1)
    # points_per_cell = typemax(T)+1 / len, but at least 1, for
    # cases where len > typemax(T)+1
    points_per_cell = one(T) << max(0, 8*sizeof(T) - bitshift)
    probs = zeros(T, at.length)
    for (i, (prob, alias)) in enumerate(at.probability_alias)
        # For intertype hasing and equality purposes we sometimes store prob = 0.5
        # even when that is above points_per_cell. Bound here:
        prob2 = min(prob, points_per_cell)
        probs[I(i-1) + alias + one(I)] += prob2
        keep = points_per_cell - prob2
        iszero(keep) || (probs[i] += keep)
        # When len > typemax(T)+1, the excess elements exist only for intertype
        # hashing and equality purposes and should be ignored.
        i > typemax(T) && break
    end
    if all(iszero, probs)
        probs[sample(zero(T), at)] = typemax(T)
    end
    probs
end

"""
    probabilities(float, at::AliasTable{T}) -> Vector{<:AbstractFloat}

Return the sampling probabilities of `at`. The returned vector will sum to 1.0, up to
rounding error.

# Example

```jldoctest
julia> AliasTables.probabilities(float, AliasTable([1, 3, 1]))
3-element Vector{Float64}:
 0.2
 0.6
 0.2
```
"""
probabilities(::typeof(float), at::AliasTable{T}) where T =
    probabilities(at) ./ (float(typemax(T))+1)

### Length accessor
"""
    length(at::AliasTable)

Get the number of weights that `at` was constructed with, including trailing zeros.

# Example

```jldoctest
julia> length(AliasTable([1, 3, 0]))
3
```
"""
Base.length(at::AliasTable) = at.length

### Show
function Base.show(io::IO, at::AliasTable{T, I}) where {T, I}
    print(io, AliasTable)
    if I != Int
        print(io, "{", T, ", ", I, "}")
    elseif T != UInt64
        print(io, "{", T, "}")
    end
    print(io, "(")
    print(IOContext(io, :typeinfo=>Vector{T}), probabilities(at))
    print(io, ")")
end
isdefined(Base, :typeinfo_implicit) && (Base.typeinfo_implicit(::Type{<:AliasTable}) = true)

### Equality and hashing

# These naive implementations are equivalent to computing equality
# based on probabilities because the constrors are deterministic w.r.t
# the input weights exclusind trailing zeros and the length is tracked.
Base.:(==)(at1::AliasTable{T, I}, at2::AliasTable{T, I}) where {T, I} =
    at1.length == at2.length && at1.probability_alias == at2.probability_alias
function Base.:(==)(at1::AliasTable{T1, I1}, at2::AliasTable{T2, I2}) where {T1, T2, I1, I2}
    at1.length == at2.length || return false
    length(at1.probability_alias) == length(at2.probability_alias) || return false
    bitshift = 8(sizeof(T1) - sizeof(T2))
    for (i,(pa1, pa2)) in enumerate(zip(at1.probability_alias, at2.probability_alias))
        pa1[2]+I1(i-1)+one(I2) == pa2[2]+I2(i-1)+one(I2) &&
        if bitshift > 0
            pa1[1] == T1(pa2[1]) << bitshift
        else
            T2(pa1[1]) << -bitshift == pa2[1]
        end || return false
    end
    true
end
struct EnumerateMapVector{T, F, P} <: AbstractVector{T}
    parent::P
    f::F
end
Base.size(mv::EnumerateMapVector) = size(mv.parent)
Base.getindex(mv::EnumerateMapVector{T, F, P}, i) where {T, F, P} = mv.f(i, mv.parent[i])
function Base.hash(at::AliasTable, h::UInt)
    h ⊻= Sys.WORD_SIZE == 32 ? 0x7719cd5e : 0x0a0c5cfeeb10f090
    h = hash(at.length, h)
    # isempty(at.probability_alias) && return hash(0, h) # This should never happen, but it makes first not throw.
    pa1 = first(at.probability_alias)
    norm(i, x) = (ldexp(float(x[1]), -8sizeof(x[1])), x[2]+typeof(x[2])(i-1)+typeof(x[2])(1))
    hash(EnumerateMapVector{typeof(norm(1, pa1)), typeof(norm), typeof(at.probability_alias)}(at.probability_alias, norm), h)
end

## Normalization

maybe_unsigned(x) = x # this is necessary to avoid calling unsigned on BigInt, Flaot64, etc.
maybe_unsigned(x::Base.BitSigned) = unsigned(x)

maybe_add_with_overflow(x::Base.BitInteger, y::Base.BitInteger) = Base.Checked.add_with_overflow(x, convert(typeof(x), y))
maybe_add_with_overflow(x, y) = x+y, false

####

# 2-4 passes (skip first two if nomralize = false)
# Initial sum, check for overflow, negatives, allzero & exactness
# If not exact, compute the mapped sum and the error there
# Dual pass for construction

# First pass
widen_to_word(x) = x
widen_to_word(x::Base.BitSignedSmall) = Int(x)
widen_to_word(x::Base.BitUnsignedSmall) = Int(x)
sum_prepare(x) = maybe_unsigned(widen_to_word(x))

function checked_sum(weights)
    xi = iterate(weights)
    xi === nothing && throw(ArgumentError("weights must be non-empty"))
    x, i = xi
    nonzero_index = 1
    while iszero(x)
        nonzero_index += 1
        xi = iterate(weights, i)
        xi === nothing && throw(ArgumentError("all weights are zero"))
        x, i = xi
    end
    x < 0 && throw(ArgumentError("found negative weight $x"))
    x0 = x
    while true
        xi = iterate(weights, i)
        xi === nothing && return (true, nonzero_index)
        x, i = xi
        iszero(x) || break
    end
    # There are two nonzero elements
    x < 0 && throw(ArgumentError("found negative weight $x"))
    sm, overflow = maybe_add_with_overflow(sum_prepare(x0), sum_prepare(x))
    overflow && !iszero(sm) && throw(ArgumentError("sum(weights) overflows"))
    while true
        xi = iterate(weights, i)
        xi === nothing && break
        x, i = xi
        x < 0 && throw(ArgumentError("found negative weight $x"))
        sm, o = maybe_add_with_overflow(sm, sum_prepare(x))
        overflow |= o
        overflow && !iszero(sm) && throw(ArgumentError("sum(weights) overflows"))
    end
    isfinite(sm) || throw(ArgumentError("sum(weights) == $sm which is not finite"))
    (false, sm)
end

# Second phase

# Slower than allocating:
# function normalize_to_uint_lazy_frac_div(::Type{T}, weights, sm) where T <: Unsigned
#     sm2 = zero(T)
#     for x in weights
#         sm2 += frac_div(T(x), sm)
#     end
#     sm2_copy = sm2 # lolz https://github.com/JuliaLang/julia/issues/15276

#     bonus = typemax(sm2)-sm2+1
#     (frac_div(T(x), sm) + (i <= bonus) for (i,x) in enumerate(weights))
# end

widen_float(T, x) = typemax(T) < floatmax(x) ? x : widen_float(T, widen(x))
function normalize_to_uint!(res::AbstractVector{T}, v, sm) where {T <: Unsigned}
    if sm isa AbstractFloat
        shift = 8sizeof(T)-exponent(sm + sqrt(eps(sm)))-1
        for i in eachindex(res, v)
            # Rounding up could cause sm2 to overflow in AliasTable(vcat(1.0-sqrt(eps(1.0)), fill(1e-100, 2^38))
            res[i] = ceil(T, ldexp(widen_float(T, v[i]), shift))
        end
        v2 = res
        sm2 = sum(res)
    else
        v2 = v
        sm2 = sm
    end

    sm3 = zero(T)

    T2 = widen(promote_type(T, typeof(sm2)))
    for (i,x) in enumerate(v2)
        # @assert x < sm2
        # @assert sm2 != 0
        val = div(T2(maybe_unsigned(x)) << 8sizeof(T), sm2) % T
        sm3 += val
        res[i] = val
    end

    sm3 == 0 && any(!iszero(res)) && return res

    for i in sm3:typemax(sm3)
        res[typemax(sm3)-i+1] += true
    end

    res
end

end
