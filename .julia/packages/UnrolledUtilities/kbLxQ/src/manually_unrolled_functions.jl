@inline _unrolled_map(::Val{0}, f, itr) = ()
@inline _unrolled_map(::Val{1}, f, itr) = (f(generic_getindex(itr, 1)),)
@inline _unrolled_map(::Val{2}, f, itr) =
    (f(generic_getindex(itr, 1)), f(generic_getindex(itr, 2)))
@inline _unrolled_map(::Val{3}, f, itr) = (
    f(generic_getindex(itr, 1)),
    f(generic_getindex(itr, 2)),
    f(generic_getindex(itr, 3)),
)
@generated _unrolled_map(::Val{N}, f, itr) where {N} =
    :(@inline Base.Cartesian.@ntuple $N n -> f(generic_getindex(itr, n)))

@inline _unrolled_any(::Val{0}, f, itr) = false
@inline _unrolled_any(::Val{1}, f, itr) = f(generic_getindex(itr, 1))
@inline _unrolled_any(::Val{2}, f, itr) =
    f(generic_getindex(itr, 1)) || f(generic_getindex(itr, 2))
@inline _unrolled_any(::Val{3}, f, itr) =
    f(generic_getindex(itr, 1)) ||
    f(generic_getindex(itr, 2)) ||
    f(generic_getindex(itr, 3))
@generated _unrolled_any(::Val{N}, f, itr) where {N} =
    :(@inline Base.Cartesian.@nany $N n -> f(generic_getindex(itr, n)))

@inline _unrolled_all(::Val{0}, f, itr) = true
@inline _unrolled_all(::Val{1}, f, itr) = f(generic_getindex(itr, 1))
@inline _unrolled_all(::Val{2}, f, itr) =
    f(generic_getindex(itr, 1)) && f(generic_getindex(itr, 2))
@inline _unrolled_all(::Val{3}, f, itr) =
    f(generic_getindex(itr, 1)) &&
    f(generic_getindex(itr, 2)) &&
    f(generic_getindex(itr, 3))
@generated _unrolled_all(::Val{N}, f, itr) where {N} =
    :(@inline Base.Cartesian.@nall $N n -> f(generic_getindex(itr, n)))

@inline _unrolled_foreach(::Val{0}, f, itr) = nothing
@inline function _unrolled_foreach(::Val{1}, f, itr)
    f(generic_getindex(itr, 1))
    return nothing
end
@inline function _unrolled_foreach(::Val{2}, f, itr)
    f(generic_getindex(itr, 1))
    f(generic_getindex(itr, 2))
    return nothing
end
@inline function _unrolled_foreach(::Val{3}, f, itr)
    f(generic_getindex(itr, 1))
    f(generic_getindex(itr, 2))
    f(generic_getindex(itr, 3))
    return nothing
end
@generated _unrolled_foreach(::Val{N}, f, itr) where {N} = quote
    @inline
    Base.Cartesian.@nexprs $N n -> f(generic_getindex(itr, n))
    return nothing
end

@inline _unrolled_reduce(::Val{0}, op, itr, init) = empty_reduction_value(init)
@inline _unrolled_reduce(::Val{1}, op, itr, init) =
    first_reduction_value(op, itr, init)
@inline _unrolled_reduce(::Val{2}, op, itr, init) =
    op(first_reduction_value(op, itr, init), generic_getindex(itr, 2))
@inline _unrolled_reduce(::Val{3}, op, itr, init) = op(
    op(first_reduction_value(op, itr, init), generic_getindex(itr, 2)),
    generic_getindex(itr, 3),
)
@generated _unrolled_reduce(::Val{N}, op, itr, init) where {N} = quote
    @inline
    value_1 = first_reduction_value(op, itr, init)
    Base.Cartesian.@nexprs $(N - 1) n ->
        (value_{n + 1} = op(value_n, generic_getindex(itr, n + 1)))
    return $(Symbol(:value_, N))
end

@inline _unrolled_accumulate(::Val{0}, op, itr, init) = ()
@inline function _unrolled_accumulate(::Val{1}, op, itr, init)
    value_1 = first_reduction_value(op, itr, init)
    return (value_1,)
end
@inline function _unrolled_accumulate(::Val{2}, op, itr, init)
    value_1 = first_reduction_value(op, itr, init)
    value_2 = op(value_1, generic_getindex(itr, 2))
    return (value_1, value_2)
end
@inline function _unrolled_accumulate(::Val{3}, op, itr, init)
    value_1 = first_reduction_value(op, itr, init)
    value_2 = op(value_1, generic_getindex(itr, 2))
    value_3 = op(value_2, generic_getindex(itr, 3))
    return (value_1, value_2, value_3)
end
@generated _unrolled_accumulate(::Val{N}, op, itr, init) where {N} = quote
    @inline
    value_1 = first_reduction_value(op, itr, init)
    Base.Cartesian.@nexprs $(N - 1) n ->
        (value_{n + 1} = op(value_n, generic_getindex(itr, n + 1)))
    return Base.Cartesian.@ntuple $N n -> value_n
end

@inline _unrolled_ifelse(::Val{0}, f, get_if, get_else, itr) = get_else()
@inline function _unrolled_ifelse(::Val{1}, f, get_if, get_else, itr)
    item_1 = generic_getindex(itr, 1)
    f(item_1) && return get_if(item_1)
    return get_else()
end
@inline function _unrolled_ifelse(::Val{2}, f, get_if, get_else, itr)
    item_1 = generic_getindex(itr, 1)
    f(item_1) && return get_if(item_1)
    item_2 = generic_getindex(itr, 2)
    f(item_2) && return get_if(item_2)
    return get_else()
end
@inline function _unrolled_ifelse(::Val{3}, f, get_if, get_else, itr)
    item_1 = generic_getindex(itr, 1)
    f(item_1) && return get_if(item_1)
    item_2 = generic_getindex(itr, 2)
    f(item_2) && return get_if(item_2)
    item_3 = generic_getindex(itr, 3)
    f(item_3) && return get_if(item_3)
    return get_else()
end
@generated _unrolled_ifelse(::Val{N}, f, get_if, get_else, itr) where {N} =
    quote
        @inline
        Base.Cartesian.@nexprs $N n -> begin
            item = generic_getindex(itr, n)
            f(item) && return get_if(item)
        end
        return get_else()
    end

@inline _unrolled_ifelse(::Val{0}, f, get_if, get_else, itr1, itr2) = get_else()
@inline function _unrolled_ifelse(::Val{1}, f, get_if, get_else, itr1, itr2)
    f(generic_getindex(itr1, 1)) && return get_if(generic_getindex(itr2, 1))
    return get_else()
end
@inline function _unrolled_ifelse(::Val{2}, f, get_if, get_else, itr1, itr2)
    f(generic_getindex(itr1, 1)) && return get_if(generic_getindex(itr2, 1))
    f(generic_getindex(itr1, 2)) && return get_if(generic_getindex(itr2, 2))
    return get_else()
end
@inline function _unrolled_ifelse(::Val{3}, f, get_if, get_else, itr1, itr2)
    f(generic_getindex(itr1, 1)) && return get_if(generic_getindex(itr2, 1))
    f(generic_getindex(itr1, 2)) && return get_if(generic_getindex(itr2, 2))
    f(generic_getindex(itr1, 3)) && return get_if(generic_getindex(itr2, 3))
    return get_else()
end
@generated _unrolled_ifelse(
    ::Val{N},
    f,
    get_if,
    get_else,
    itr1,
    itr2,
) where {N} = quote
    @inline
    Base.Cartesian.@nexprs $N n -> begin
        f(generic_getindex(itr1, n)) && return get_if(generic_getindex(itr2, n))
    end
    return get_else()
end
