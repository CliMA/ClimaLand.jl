"""
    StaticBitVector{N, [U]}(f)
    StaticBitVector{N, [U]}([bit])

A statically sized analogue of `BitVector` with `Unsigned` chunks of type `U`,
which can be constructed using either a function `f(n)` or a constant `bit`. By
default, `U` is set to `UInt8` and `bit` is set to `false`.

This iterator can only store `Bool`s, so its `output_type_for_promotion` is a
`ConditionalOutputType`. Efficient implementations are provided for all unrolled
functions, though the methods for `unrolled_map` and `unrolled_accumulate` only
apply when the first item in the output is a `Bool`.
"""
struct StaticBitVector{N, U <: Unsigned, I <: NTuple{<:Any, U}} <:
       StaticSequence{N}
    ints::I
end
@inline StaticBitVector{N, U}(ints) where {N, U} =
    StaticBitVector{N, U, typeof(ints)}(ints)
@inline StaticBitVector{N}(args...) where {N} =
    StaticBitVector{N, UInt8}(args...)

@inline function StaticBitVector{N, U}(bit::Bool = false) where {N, U}
    n_bits_per_int = 8 * sizeof(U)
    n_ints = cld(N, n_bits_per_int)
    ints = ntuple(Returns(bit ? ~zero(U) : zero(U)), Val(n_ints))
    return StaticBitVector{N, U}(ints)
end

@inline function StaticBitVector{N, U}(f::Function) where {N, U}
    n_bits_per_int = 8 * sizeof(U)
    n_ints = cld(N, n_bits_per_int)
    ints = ntuple(Val(n_ints)) do int_index
        @inline
        first_index = n_bits_per_int * (int_index - 1) + 1
        unrolled_reduce(
            StaticOneTo(min(n_bits_per_int, N - first_index + 1)),
            zero(U),
        ) do int, bit_index
            @inline
            bit_offset = bit_index - 1
            int | U(f(first_index + bit_offset)::Bool) << bit_offset
        end
    end
    return StaticBitVector{N, U}(ints)
end

@inline function int_index_and_bit_offset(::Type{U}, n) where {U}
    int_offset, bit_offset = divrem(n - 1, 8 * sizeof(U))
    return (int_offset + 1, bit_offset)
end

@inline function generic_getindex(
    itr::StaticBitVector{<:Any, U},
    n::Integer,
) where {U}
    int_index, bit_offset = int_index_and_bit_offset(U, n)
    int = itr.ints[int_index]
    return Bool(int >> bit_offset & one(int))
end

@inline function Base.setindex(
    itr::StaticBitVector{N, U},
    bit::Bool,
    n::Integer,
) where {N, U}
    int_index, bit_offset = int_index_and_bit_offset(U, n)
    int = itr.ints[int_index]
    new_int = int & ~(one(U) << bit_offset) | U(bit) << bit_offset
    ints = Base.setindex(itr.ints, new_int, int_index)
    return StaticBitVector{N, U}(ints)
end

@inline output_type_for_promotion(::StaticBitVector{<:Any, U}) where {U} =
    ConditionalOutputType(Bool, StaticBitVector{<:Any, U})

@inline empty_output(::Type{StaticBitVector{<:Any, U}}) where {U} =
    StaticBitVector{0, U}()

@inline unrolled_map_into(::Type{StaticBitVector{<:Any, U}}, f, itr) where {U} =
    StaticBitVector{length(itr), U}(
        Base.Fix1(generic_getindex, Iterators.map(f, itr)),
    )

@inline function unrolled_push_into(
    ::Type{StaticBitVector{<:Any, U}},
    itr,
    bit,
) where {U}
    n_bits_per_int = 8 * sizeof(U)
    n_ints = cld(length(itr), n_bits_per_int)
    bit_offset = length(itr) % n_bits_per_int
    ints = if bit_offset == 0
        (itr.ints..., U(bit))
    else
        last_int = itr.ints[n_ints]
        new_last_int =
            last_int & ~(one(U) << bit_offset) | U(bit) << bit_offset
        (unrolled_take(itr.ints, Val(n_ints - 1))..., new_last_int)
    end
    return StaticBitVector{length(itr) + 1, U}(ints)
end

@inline function unrolled_append_into(
    ::Type{StaticBitVector{<:Any, U}},
    itr1,
    itr2,
) where {U}
    n_bits_per_int = 8 * sizeof(U)
    n_ints1 = cld(length(itr1), n_bits_per_int)
    bit_offset = length(itr1) % n_bits_per_int
    ints = if bit_offset == 0 || length(itr2) == 0
        (itr1.ints..., itr2.ints...)
    else
        mid_int1 = itr1.ints[n_ints1]
        mid_int2 = itr2.ints[1]
        mid_int =
            mid_int1 & ~(~zero(U) << bit_offset) | mid_int2 << bit_offset
        final_ints =
            length(itr2) + bit_offset <= n_bits_per_int ? () :
            unrolled_drop(itr2, Val(n_bits_per_int - bit_offset)).ints
        (unrolled_take(itr1.ints, Val(n_ints1 - 1))..., mid_int, final_ints...)
    end
    return StaticBitVector{length(itr1) + length(itr2), U}(ints)
end

@inline function unrolled_take_into(
    ::Type{StaticBitVector{<:Any, U}},
    itr,
    ::Val{N},
) where {N, U}
    n_bits_per_int = 8 * sizeof(U)
    n_ints = cld(N, n_bits_per_int)
    ints = unrolled_take(itr.ints, Val(n_ints))
    return StaticBitVector{N, U}(ints)
end

@inline function unrolled_drop_into(
    ::Type{StaticBitVector{<:Any, U}},
    itr,
    ::Val{N},
) where {N, U}
    n_bits_per_int = 8 * sizeof(U)
    n_ints = cld(length(itr) - N, n_bits_per_int)
    n_dropped_ints = fld(N, n_bits_per_int)
    bit_offset = N - n_bits_per_int * n_dropped_ints
    ints_without_offset = unrolled_drop(itr.ints, Val(n_dropped_ints))
    ints = if bit_offset == 0 || length(itr) <= N
        ints_without_offset
    else
        next_ints =
            length(ints_without_offset) == 1 ? (nothing,) :
            (unrolled_drop(ints_without_offset, Val(1))..., nothing)
        unrolled_map(ints_without_offset, next_ints) do cur_int, next_int
            @inline
            isnothing(next_int) ? cur_int >> bit_offset :
            cur_int >> bit_offset | next_int << (n_bits_per_int - bit_offset)
        end
    end
    return StaticBitVector{length(itr) - N, U}(ints)
end

@inline function unrolled_accumulate_into(
    ::Type{StaticBitVector{<:Any, U}},
    op,
    itr,
    init,
) where {U}
    N = length(itr)
    n_bits_per_int = 8 * sizeof(U)
    n_ints = cld(N, n_bits_per_int)
    ints = unrolled_accumulate(
        StaticOneTo(n_ints),
        (nothing, init),
        first,
    ) do (_, init_value_for_new_int), int_index
        @inline
        first_index = n_bits_per_int * (int_index - 1) + 1
        unrolled_reduce(
            StaticOneTo(min(n_bits_per_int, N - first_index + 1)),
            (zero(U), init_value_for_new_int),
        ) do (int, prev_value), bit_index
            @inline
            bit_offset = bit_index - 1
            item = generic_getindex(itr, first_index + bit_offset)
            new_value =
                first_index + bit_offset == 1 && prev_value isa NoInit ?
                item : op(prev_value, item)
            (int | U(new_value::Bool) << bit_offset, new_value)
        end
    end
    return StaticBitVector{N, U}(ints)
end
