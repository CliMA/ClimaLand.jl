"""
    generic_getindex(itr, n)

Identical to `getindex(itr, n)`, but with the added ability to handle lazy
iterator types defined in the standard library, such as `Base.Generator` and
`Iterators.Enumerate`.
"""
@inline generic_getindex(itr, n) = getindex(itr, n)
@inline generic_getindex(itr::Base.Generator, n) =
    itr.f(generic_getindex(itr.iter, n))
@inline generic_getindex(itr::Iterators.Reverse, n) =
    generic_getindex(itr.itr, length(itr.itr) - n + 1)
@inline generic_getindex(itr::Iterators.Enumerate, n) =
    (n, generic_getindex(itr.itr, n))
@inline generic_getindex(itr::Iterators.Zip, n) =
    unrolled_map(Base.Fix2(generic_getindex, n), itr.is)

@inline first_item_type(itr) =
    Base.promote_op(Base.Fix2(generic_getindex, 1), typeof(itr))
@inline second_item_type(itr) =
    Base.promote_op(Base.Fix2(generic_getindex, 2), typeof(itr))

"""
    output_type_for_promotion(itr)

The type of output that unrolled functions should try to generate for the input
iterator `itr`, or a `ConditionalOutputType` if the output type depends on the
type of items that need to be stored in it, or `NoOutputType()` if `itr` is a
lazy iterator without any associated output type. Defaults to `Tuple`.
"""
@inline output_type_for_promotion(_) = Tuple
@inline output_type_for_promotion(::NamedTuple{names}) where {names} =
    NamedTuple{names}
@inline output_type_for_promotion(itr::Base.Generator) =
    output_type_for_promotion(itr.iter)
@inline output_type_for_promotion(itr::Iterators.Reverse) =
    output_type_for_promotion(itr.itr)
@inline output_type_for_promotion(itr::Iterators.Enumerate) =
    output_type_for_promotion(itr.itr)
@inline output_type_for_promotion(itr::Iterators.Zip) =
    maybe_ambiguous_promoted_output_type(itr.is...)

"""
    AmbiguousOutputType

The result of `output_type_for_promotion` for iterators that do not have
well-defined output types.
"""
abstract type AmbiguousOutputType end

"""
    NoOutputType()

The `AmbiguousOutputType` of lazy iterators.
"""
struct NoOutputType <: AmbiguousOutputType end

"""
    ConditionalOutputType(allowed_item_type, output_type, [fallback_type])

An `AmbiguousOutputType` that can have one of two possible values. If the first
item in the output is a subtype of `allowed_item_type`, the output will have the
type `output_type`; otherwise, it will have the type `fallback_type`, which is
set to `Tuple` by default.
"""
struct ConditionalOutputType{I, O, O′} <: AmbiguousOutputType end
@inline ConditionalOutputType(
    allowed_item_type::Type,
    output_type::Type,
    fallback_type::Type = Tuple,
) = ConditionalOutputType{allowed_item_type, output_type, fallback_type}()

@inline unambiguous_output_type(_, ::Type{O}) where {O} = O
@inline unambiguous_output_type(_, ::NoOutputType) = Tuple
@inline unambiguous_output_type(
    get_first_item_type,
    ::ConditionalOutputType{I, O, O′},
) where {I, O, O′} = get_first_item_type() <: I ? O : O′

"""
    output_promote_rule(output_type1, output_type2)

The type of output that should be generated when two iterators do not have the
same `output_type_for_promotion`, or `Union{}` if these iterators should not be
used together. Only one method of `output_promote_rule` needs to be defined for
any pair of output types.

By default, all types take precedence over `NoOutputType()`, and the conditional
part of any `ConditionalOutputType` takes precedence over an unconditional type
(so that only the `fallback_type` of any conditional type gets promoted). The
default result for all other pairs of unequal output types is `Union{}`.
"""
@inline output_promote_rule(_, _) = Union{}
@inline output_promote_rule(::Type{O}, ::Type{O}) where {O} = O
@inline output_promote_rule(::NoOutputType, output_type) = output_type
@inline output_promote_rule(
    ::ConditionalOutputType{I, O, O′},
    ::Type{O′′},
) where {I, O, O′, O′′} =
    ConditionalOutputType(I, O, output_promote_rule(O′, O′′))
@inline output_promote_rule(
    ::Type{O′},
    ::ConditionalOutputType{I, O, O′′},
) where {I, O, O′, O′′} =
    ConditionalOutputType(I, O, output_promote_rule(O′, O′′))
@inline output_promote_rule(
    ::ConditionalOutputType{I, O, O′},
    ::ConditionalOutputType{I, O, O′′},
) where {I, O, O′, O′′} =
    ConditionalOutputType(I, O, output_promote_rule(O′, O′′))

@inline function output_promote_result(O1, O2)
    O12 = output_promote_rule(O1, O2)
    O21 = output_promote_rule(O2, O1)
    O12 == O21 == Union{} &&
        error("output_promote_rule is undefined for $O1 and $O2")
    (O12 == O21 || O21 == Union{}) && return O12
    O12 == Union{} && return O21
    error("output_promote_rule yields inconsistent results for $O1 and $O2: \
           $O12 for $O1 followed by $O2, versus $O21 for $O2 followed by $O1")
end

@inline maybe_ambiguous_promoted_output_type(itrs...) =
    isempty(itrs) ? Tuple : # Generate a Tuple when given 0 inputs.
    unrolled_mapreduce(output_type_for_promotion, output_promote_result, itrs)

@inline inferred_output_type(itr) =
    unambiguous_output_type(output_type_for_promotion(itr)) do
        @inline
        first_item_type(itr)
    end

@inline promoted_output_type(itrs...) =
    unambiguous_output_type(maybe_ambiguous_promoted_output_type(itrs...)) do
        @inline
        first_item_type(generic_getindex(itrs, 1))
    end

@inline unrolled_map_output_type(f, itr) =
    inferred_output_type(Iterators.map(f, itr))

@inline unrolled_accumulate_output_type(op, itr, init) =
    unambiguous_output_type(output_type_for_promotion(itr)) do
        @inline
        no_init = init isa NoInit
        arg1_type = no_init ? first_item_type(itr) : typeof(init)
        arg2_type = no_init ? second_item_type(itr) : first_item_type(itr)
        Base.promote_op(op, arg1_type, arg2_type)
    end

"""
    constructor_from_tuple(output_type)

A function that can be used to efficiently construct an output of type
`output_type` from a `Tuple`, or `identity` if such an output should not be
constructed from a `Tuple`. Defaults to `identity`, which also handles the case
where `output_type` is already `Tuple`. The `output_type` here is guaranteed to
be a `Type`, rather than a `ConditionalOutputType` or `NoOutputType`.

Many statically sized iterators (e.g., `SVector`s) are essentially wrappers for
`Tuple`s, and their constructors for `Tuple`s can be reduced to no-ops. The main
exceptions are [`StaticOneTo`](@ref UnrolledUtilities.StaticOneTo)s and
[`StaticBitVector`](@ref UnrolledUtilities.StaticBitVector)s, which do not
provide constructors for `Tuple`s because there is no performance benefit to
making a lazy or low-storage data structure once a corresponding high-storage
data structure has already been constructed.
"""
@inline constructor_from_tuple(::Type) = identity
@inline constructor_from_tuple(::Type{NT}) where {NT <: NamedTuple} = NT

"""
    empty_output(output_type)

An empty output of type `output_type`. Defaults to applying the
`constructor_from_tuple` for the given type to an empty `Tuple`.
"""
@inline empty_output(output_type) = constructor_from_tuple(output_type)(())

@inline inferred_empty(itr) = empty_output(inferred_output_type(itr))

# This makes lazy iterators non-lazy, and it is a no-op for non-lazy iterators.
@inline non_lazy_iterator(itr) = unrolled_append(itr, inferred_empty(itr))
