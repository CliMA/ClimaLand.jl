module UnrolledUtilities

export StaticSequence,
    StaticOneTo,
    StaticBitVector,
    unrolled_push,
    unrolled_append,
    unrolled_prepend,
    unrolled_take,
    unrolled_drop,
    unrolled_map,
    unrolled_any,
    unrolled_all,
    unrolled_foreach,
    unrolled_reduce,
    unrolled_mapreduce,
    unrolled_accumulate,
    unrolled_applyat,
    unrolled_in,
    unrolled_unique,
    unrolled_allunique,
    unrolled_allequal,
    unrolled_sum,
    unrolled_prod,
    unrolled_cumsum,
    unrolled_cumprod,
    unrolled_count,
    unrolled_maximum,
    unrolled_minimum,
    unrolled_extrema,
    unrolled_findmax,
    unrolled_findmin,
    unrolled_argmax,
    unrolled_argmin,
    unrolled_findfirst,
    unrolled_findlast,
    unrolled_argfirst,
    unrolled_arglast,
    unrolled_filter,
    unrolled_split,
    unrolled_flatten,
    unrolled_flatmap,
    unrolled_product,
    unrolled_cycle,
    unrolled_partition

##
## Internal types and functions
##

include("unrollable_iterator_interface.jl")

# Analogue of the non-public Base._InitialValue for reduction and accumulation.
struct NoInit end

@inline empty_reduction_value(::NoInit) =
    error("unrolled_reduce requires an init value for empty iterators")
@inline empty_reduction_value(init) = init

@inline first_reduction_value(op, itr, ::NoInit) = generic_getindex(itr, 1)
@inline first_reduction_value(op, itr, init) =
    op(init, generic_getindex(itr, 1))

# Analogue of ∘, but with only one function argument and guaranteed inlining.
# Base's ∘ leads to type instabilities in unit tests on Julia 1.10 and 1.11.
@inline ⋅(f1::F1, f2::F2) where {F1, F2} = x -> (@inline f1(f2(x)))
@inline ⋅(::Type{T}, f::F) where {T, F} = x -> (@inline T(f(x)))
@inline ⋅(f::F, ::Type{T}) where {F, T} = x -> (@inline f(T(x)))
@inline ⋅(::Type{T1}, ::Type{T2}) where {T1, T2} = x -> (@inline T1(T2(x)))

##
## Types with optimized methods for unrolled functions
##

"""
    StaticSequence{N}

Abstract type that can represent any iterable of constant length `N`. Subtypes
include the low-storage data structures `StaticOneTo` and `StaticBitVector`,
which have optimized methods for certain unrolled functions.
"""
abstract type StaticSequence{N} end

@inline Base.length(::StaticSequence{N}) where {N} = N
@inline Base.firstindex(::StaticSequence) = 1
@inline Base.lastindex(itr::StaticSequence) = length(itr)
@inline Base.getindex(itr::StaticSequence, n::Integer) =
    generic_getindex(itr, n)
@inline Base.iterate(itr::StaticSequence, n = 1) =
    n > length(itr) ? nothing : (generic_getindex(itr, n), n + 1)
@inline Base.iterate(
    itr::Iterators.Reverse{<:StaticSequence},
    n = length(itr),
) = n < 1 ? nothing : (generic_getindex(itr.itr, n), n - 1)

include("StaticOneTo.jl")
include("StaticBitVector.jl")

@inline static_range(itr) = StaticOneTo(length(itr))

##
## Functions unrolled using either Core.tuple or ntuple
##

@inline unrolled_push_into(output_type, itr, item) =
    constructor_from_tuple(output_type)((itr..., item))
@inline unrolled_push(itr, item) =
    unrolled_push_into(inferred_output_type(itr), itr, item)
@inline unrolled_push(itr, items...) =
    unrolled_reduce(unrolled_push, items, itr)

@inline unrolled_append_into(output_type, itr1, itr2) =
    constructor_from_tuple(output_type)(
        ntuple(Val(length(itr1) + length(itr2))) do n
            @inline
            n <= length(itr1) ? generic_getindex(itr1, n) :
            generic_getindex(itr2, n - length(itr1))
        end,
    )
@inline unrolled_append(itr1, itr2) =
    unrolled_append_into(promoted_output_type(itr1, itr2), itr1, itr2)
@inline unrolled_append(itr, itrs...) =
    unrolled_reduce(unrolled_append, itrs, itr)

@inline unrolled_prepend(itr, itrs...) = unrolled_append(itrs..., itr)

@inline unrolled_take_into(output_type, itr, ::Val{N}) where {N} =
    constructor_from_tuple(output_type)(
        ntuple(Base.Fix1(generic_getindex, itr), Val(N)),
    )
@inline unrolled_take(itr, val_N) =
    unrolled_take_into(inferred_output_type(itr), itr, val_N)

@inline unrolled_drop_into(output_type, itr, ::Val{N}) where {N} =
    constructor_from_tuple(output_type)(
        ntuple(
            Base.Fix1(generic_getindex, itr) ⋅ Base.Fix1(+, N),
            Val(length(itr) - N),
        ),
    )
@inline unrolled_drop(itr, val_N) =
    unrolled_drop_into(inferred_output_type(itr), itr, val_N)

##
## Functions unrolled using either hard-coded or generated expressions
##

include("manually_unrolled_functions.jl")

# The unrolled_map function could also be implemented in terms of ntuple, but
# then it would be subject to the same recursion limit as ntuple. On Julia 1.10,
# this leads to type instabilities in several unit tests for nested iterators.
@inline unrolled_map_into_tuple(f::F, itr) where {F} =
    _unrolled_map(Val(length(itr)), f, itr)
@inline unrolled_map_into(output_type, f::F, itr) where {F} =
    constructor_from_tuple(output_type)(unrolled_map_into_tuple(f, itr))
@inline unrolled_map(f::F, itr) where {F} =
    unrolled_map_into(inferred_output_type(Iterators.map(f, itr)), f, itr)
@inline unrolled_map(f::F, itrs...) where {F} =
    unrolled_map(splat(f), zip(itrs...))

@inline unrolled_any(itr) = unrolled_any(identity, itr)
@inline unrolled_any(f::F, itr) where {F} =
    _unrolled_any(Val(length(itr)), f, itr)

@inline unrolled_all(itr) = unrolled_all(identity, itr)
@inline unrolled_all(f::F, itr) where {F} =
    _unrolled_all(Val(length(itr)), f, itr)

@inline unrolled_foreach(f::F, itr) where {F} =
    _unrolled_foreach(Val(length(itr)), f, itr)
@inline unrolled_foreach(f, itrs...) = unrolled_foreach(splat(f), zip(itrs...))

@inline unrolled_reduce(op::O, itr, init) where {O} =
    _unrolled_reduce(Val(length(itr)), op, itr, init)
@inline unrolled_reduce(op::O, itr; init = NoInit()) where {O} =
    unrolled_reduce(op, itr, init)

@inline unrolled_mapreduce(f::F, op::O, itrs...; init = NoInit()) where {F, O} =
    unrolled_reduce(op, unrolled_map(f, itrs...), init)

@inline unrolled_accumulate_into_tuple(op::O, itr, init) where {O} =
    _unrolled_accumulate(Val(length(itr)), op, itr, init)
@inline unrolled_accumulate_into(output_type, op::O, itr, init) where {O} =
    constructor_from_tuple(output_type)(
        unrolled_accumulate_into_tuple(op, itr, init),
    )
@inline unrolled_accumulate(op::O, itr, init) where {O} =
    unrolled_accumulate_into(
        unrolled_accumulate_output_type(op, itr, init),
        op,
        itr,
        init,
    )
@inline unrolled_accumulate(op::O, itr; init = NoInit()) where {O} =
    unrolled_accumulate(op, itr, init)

# The unrolled_ifelse function is for internal use, and it is not exported.
# With one iterator as an argument, it is similar to
# f(itr[1]) ? get_if(itr[1]) : f(itr[2]) ? get_if(itr[2]) ... : get_else().
# With two iterators as arguments, it is similar to
# f(itr1[1]) ? get_if(itr2[1]) : f(itr1[2]) ? get_if(itr2[2]) ... : get_else().
# When f compares a constant isbits value computed from each item against a
# non-constant isbits value (as in unrolled_applyat), unrolled_ifelse is
# optimized into a switch instruction during LLVM code generation.
@inline unrolled_ifelse(
    f::F,
    get_if::I,
    get_else::E,
    itr,
    itrs...,
) where {F, I, E} =
    _unrolled_ifelse(Val(length(itr)), f, get_if, get_else, itr, itrs...)

##
## Unrolled functions without any analogues in Base
##

@inline unrolled_applyat(f::F, n, itr) where {F} = unrolled_ifelse(
    ==(n),
    f,
    () -> throw(BoundsError(itr, n)),
    static_range(itr),
    itr,
)
@inline unrolled_applyat(f::F, n, itrs...) where {F} =
    unrolled_applyat(splat(f), n, zip(itrs...))

##
## Unrolled analogues of functions from base/operators.jl and base/set.jl
##

# Using === instead of == or isequal improves type stability for singletons.

@inline unrolled_in(item, itr) = unrolled_any(Base.Fix1(===, item), itr)

@inline unrolled_unique(itr) =
    unrolled_reduce(itr, inferred_empty(itr)) do items, item
        @inline
        unrolled_in(item, items) ? items : unrolled_push(items, item)
    end
@inline unrolled_unique(f::F, itr) where {F} =
    unrolled_reduce(itr, ((), inferred_empty(itr))) do (f_values, items), item
        @inline
        f_value = f(item)
        unrolled_in(f_value, f_values) ? (f_values, items) :
        (unrolled_push(f_values, f_value), unrolled_push(items, item))
    end[2]

@inline unrolled_allunique(itr) = unrolled_allunique(identity, itr)
@inline unrolled_allunique(f::F, itr) where {F} =
    length(itr) == length(unrolled_unique(f, itr))

@inline unrolled_allequal(itr) = unrolled_allequal(identity, itr)
@inline unrolled_allequal(f::F, itr) where {F} =
    isempty(itr) ? true :
    unrolled_all(
        Base.Fix1(===, f(generic_getindex(itr, 1))) ⋅ f,
        unrolled_drop(itr, Val(1)),
    )

##
## Unrolled analogues of functions from base/reduce.jl and base/accumulate.jl
##

# Sum and prod only need init when itr is empty, so it can be ignored otherwise.

@inline unrolled_sum(itr; init = 0) = unrolled_sum(identity, itr; init)
@inline unrolled_sum(f::F, itr; init = 0) where {F} =
    isempty(itr) ? init : unrolled_mapreduce(f, +, itr)

@inline unrolled_prod(itr; init = 1) = unrolled_prod(identity, itr; init)
@inline unrolled_prod(f::F, itr; init = 1) where {F} =
    isempty(itr) ? init : unrolled_mapreduce(f, *, itr)

@inline unrolled_cumsum(itr) = unrolled_cumsum(identity, itr)
@inline unrolled_cumsum(f::F, itr) where {F} =
    unrolled_accumulate(+, unrolled_map(f, itr))

@inline unrolled_cumprod(itr) = unrolled_cumprod(identity, itr)
@inline unrolled_cumprod(f::F, itr) where {F} =
    unrolled_accumulate(*, unrolled_map(f, itr))

@inline unrolled_count(itr) = unrolled_count(identity, itr)
@inline unrolled_count(f::F, itr) where {F} = unrolled_sum(Bool ⋅ f, itr)

@inline unrolled_maximum(itr) = unrolled_maximum(identity, itr)
@inline unrolled_maximum(f::F, itr) where {F} = unrolled_mapreduce(f, max, itr)

@inline unrolled_minimum(itr) = unrolled_minimum(identity, itr)
@inline unrolled_minimum(f::F, itr) where {F} = unrolled_mapreduce(f, min, itr)

@inline extrema_reduction_operator((f_min1, f_max1), (f_min2, f_max2)) =
    (min(f_min1, f_min2), max(f_max1, f_max2))
@inline unrolled_extrema(itr) = unrolled_extrema(identity, itr)
@inline unrolled_extrema(f::F, itr) where {F} =
    unrolled_mapreduce(extrema_reduction_operator, itr) do item
        @inline
        f_value = f(item)
        (f_value, f_value)
    end

@inline findmax_reduction_operator((f_value1, value1), (f_value2, value2)) =
    f_value1 < f_value2 ? (f_value2, value2) : (f_value1, value1)
@inline unrolled_findmax(itr) = unrolled_findmax(identity, itr)
@inline unrolled_findmax(f::F, itr) where {F} =
    unrolled_mapreduce(findmax_reduction_operator, enumerate(itr)) do (n, item)
        @inline
        (f(item), n)
    end

@inline findmin_reduction_operator((f_value1, value1), (f_value2, value2)) =
    f_value1 > f_value2 ? (f_value2, value2) : (f_value1, value1)
@inline unrolled_findmin(itr) = unrolled_findmin(identity, itr)
@inline unrolled_findmin(f::F, itr) where {F} =
    unrolled_mapreduce(findmin_reduction_operator, enumerate(itr)) do (n, item)
        @inline
        (f(item), n)
    end

@inline unrolled_argmax(itr) = unrolled_findmax(itr)[2]
@inline unrolled_argmax(f::F, itr) where {F} =
    unrolled_mapreduce(findmax_reduction_operator, itr) do item
        @inline
        (f(item), item)
    end[2]

@inline unrolled_argmin(itr) = unrolled_findmin(itr)[2]
@inline unrolled_argmin(f::F, itr) where {F} =
    unrolled_mapreduce(findmin_reduction_operator, itr) do item
        @inline
        (f(item), item)
    end[2]

##
## Unrolled analogues of functions from base/arrays.jl
##

@inline unrolled_findfirst(itr) = unrolled_findfirst(identity, itr)
@inline unrolled_findfirst(f::F, itr) where {F} =
    unrolled_ifelse(f, identity, Returns(nothing), itr, static_range(itr))

@inline unrolled_findlast(itr) = unrolled_findlast(identity, itr)
@inline unrolled_findlast(f::F, itr) where {F} = unrolled_ifelse(
    f,
    identity,
    Returns(nothing),
    Iterators.reverse(itr),
    Iterators.reverse(static_range(itr)),
)

@inline unrolled_argfirst(f::F, itr) where {F} = unrolled_ifelse(
    f,
    identity,
    () -> error("itr does not contain any items for which f(item) is true"),
    itr,
)

@inline unrolled_arglast(f::F, itr) where {F} =
    unrolled_argfirst(f, Iterators.reverse(itr))

@inline unrolled_filter(f::F, itr) where {F} =
    unrolled_reduce(itr, inferred_empty(itr)) do items_with_true_f, item
        @inline
        f(item) ? unrolled_push(items_with_true_f, item) : items_with_true_f
    end

@inline unrolled_split(f::F, itr) where {F} =
    unrolled_reduce(
        itr,
        (inferred_empty(itr), inferred_empty(itr)),
    ) do (items_with_true_f, items_with_false_f), item
        @inline
        f(item) ? (unrolled_push(items_with_true_f, item), items_with_false_f) :
        (items_with_true_f, unrolled_push(items_with_false_f, item))
    end

##
## Unrolled analogues of functions from base/iterators.jl
##

@inline unrolled_flatten(itr) =
    if isempty(itr)
        inferred_empty(itr)
    elseif length(itr) == 1
        non_lazy_iterator(generic_getindex(itr, 1))
    else
        unrolled_reduce(unrolled_append, itr)
    end

@inline unrolled_flatmap(f::F, itrs...) where {F} =
    unrolled_flatten(unrolled_map(f, itrs...))

@inline unrolled_product(itrs...) =
    ntuple(Val(unrolled_prod(length, itrs))) do n
        @inline
        Base.@assume_effects :foldable
        items = ntuple(Val(length(itrs))) do itr_index
            @inline
            Base.@assume_effects :foldable
            cur_length = length(itrs[itr_index])
            prev_length = unrolled_prod(length, itrs[1:(itr_index - 1)])
            item_index = (n - 1) ÷ prev_length % cur_length + 1
            generic_getindex(itrs[itr_index], item_index)
        end
        # Since the iterators are passed as individual arguments, the number of
        # items is limited by the number of argument registers, so we can just
        # use constructor_from_tuple and avoid sacrificing compilation time to
        # optimize for low-storage iterators (e.g., by using unrolled_push).
        constructor_from_tuple(promoted_output_type(itrs...))(items)
    end

@inline unrolled_cycle(itr, ::Val{N}) where {N} =
    unrolled_flatten(ntuple(Returns(itr), Val(N)))

@inline unrolled_partition(itr, ::Val{N}) where {N} =
    ntuple(Val(cld(length(itr), N))) do partition_number
        @inline
        first_index = N * (partition_number - 1)
        last_index = min(length(itr), N * partition_number)
        unrolled_drop(unrolled_take(itr, Val(last_index)), Val(first_index))
    end

##
## Additional instructions for compilation
##

# Remove the default recursion limit from every function defined in this module.
@static if hasfield(Method, :recursion_relation)
    module_names = names(@__MODULE__; all = true)
    module_values = map(Base.Fix1(getproperty, @__MODULE__), module_names)
    module_functions = filter(Base.Fix2(isa, Function), module_values)
    for f in module_functions, method in methods(f)
        method.recursion_relation = Returns(true)
    end
end

end
