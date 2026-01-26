using UnrolledUtilities: NoInit, generic_getindex, unrolled_drop

@inline _rec_unrolled_map(f) = ()
@inline _rec_unrolled_map(f, item, items...) =
    (f(item), _rec_unrolled_map(f, items...)...)
@inline rec_unrolled_map(f, itr) = _rec_unrolled_map(f, itr...)

@inline _rec_unrolled_any(f) = false
@inline _rec_unrolled_any(f, item, items...) =
    f(item) || _rec_unrolled_any(f, items...)
@inline rec_unrolled_any(f, itr) = _rec_unrolled_any(f, itr...)

@inline _rec_unrolled_all(f) = true
@inline _rec_unrolled_all(f, item, items...) =
    f(item) && _rec_unrolled_all(f, items...)
@inline rec_unrolled_all(f, itr) = _rec_unrolled_all(f, itr...)

@inline _rec_unrolled_foreach(f) = nothing
@inline _rec_unrolled_foreach(f, item, items...) =
    (f(item); _rec_unrolled_foreach(f, items...))
@inline rec_unrolled_foreach(f, itr) = _rec_unrolled_foreach(f, itr...)

@inline _rec_unrolled_reduce(op, prev_value) = prev_value
@inline _rec_unrolled_reduce(op, prev_value, item, items...) =
    _rec_unrolled_reduce(op, op(prev_value, item), items...)
@inline rec_unrolled_reduce(op, itr, init) =
    init isa NoInit ? _rec_unrolled_reduce(op, itr...) :
    _rec_unrolled_reduce(op, init, itr...)

@inline _rec_unrolled_accumulate(op, prev_value) = (prev_value,)
@inline _rec_unrolled_accumulate(op, prev_value, item, items...) = (
    prev_value,
    _rec_unrolled_accumulate(op, op(prev_value, item), items...)...,
)
@inline rec_unrolled_accumulate(op, itr, init) =
    isempty(itr) ? () :
    init isa NoInit ? _rec_unrolled_accumulate(op, itr...) :
    _rec_unrolled_accumulate(
        op,
        op(init, generic_getindex(itr, 1)),
        unrolled_drop(itr, Val(1))...,
    )

@inline _rec_unrolled_ifelse(f, get_if, get_else) = get_else()
@inline _rec_unrolled_ifelse(f, get_if, get_else, item, items...) =
    f(item) ? get_if(item) : _rec_unrolled_ifelse(f, get_if, get_else, items...)
@inline rec_unrolled_ifelse(f, get_if, get_else, itr) =
    _rec_unrolled_ifelse(f, get_if, get_else, itr...)

@inline _rec_unrolled_ifelse2(f, get_if, get_else) = get_else()
@inline _rec_unrolled_ifelse2(f, get_if, get_else, (item1, item2), items...) =
    f(item1) ? get_if(item2) :
    _rec_unrolled_ifelse2(f, get_if, get_else, items...)
@inline rec_unrolled_ifelse(f, get_if, get_else, itr1, itr2) =
    _rec_unrolled_ifelse2(f, get_if, get_else, zip(itr1, itr2)...)
