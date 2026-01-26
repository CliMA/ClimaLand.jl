```@meta
CurrentModule = UnrolledUtilities
```

## How to Unroll

There are two general ways to implement loop unrolling in Julia—recursively
splatting iterator contents and manually constructing unrolled expressions. For
example, a recursively unrolled version of the `foreach` function is

```julia
unrolled_foreach(f, itr) = _unrolled_foreach(f, itr...)
_unrolled_foreach(f) = nothing
_unrolled_foreach(f, item, items...) = (f(item); _unrolled_foreach(f, items...))
```

In contrast, a manually unrolled implementation of this function looks like

```julia
unrolled_foreach(f, itr) = _unrolled_foreach(Val(length(itr)), f, itr)
@generated _unrolled_foreach(::Val{N}, f, itr) where {N} =
    Expr(:block, (:(f(generic_getindex(itr, $n))) for n in 1:N)..., nothing)
```

Julia's compiler can only pass up to 32 values through function arguments
without allocating heap memory, so recursive unrolling is not type-stable for
iterators with lengths greater than 32. However, automatically generating
functions often requires more time and memory resources during compilation than
writing hard-coded functions. Recursive inlining adds overhead to compilation
as well, but this is typically smaller than the overhead of generated functions
for short iterators. To avoid sacrificing latency by using generated functions,
several hard-coded methods can be added to the manually unrolled implementation:

```julia
_unrolled_foreach(::Val{0}, f, itr) = nothing
_unrolled_foreach(::Val{1}, f, itr) = (f(generic_getindex(itr, 1)); nothing)
_unrolled_foreach(::Val{2}, f, itr) =
    (f(generic_getindex(itr, 1)); f(generic_getindex(itr, 2)); nothing)
_unrolled_foreach(::Val{3}, f, itr) =
    (f(generic_getindex(itr, 1)); f(generic_getindex(itr, 2)); f(generic_getindex(itr, 3)); nothing)
```

With this modification, manual unrolling does not exceed the compilation
requirements of recursive unrolling across a wide range of use cases. Since it
also avoids type instabilities for arbitrarily large iterators, a combination
of hard-coded and generated functions with manual unrolling serves as the basis
of all unrolled functions defined in this package. Similarly, the
[`ntuple(f, ::Val{N})`](https://github.com/JuliaLang/julia/blob/v1.11.0/base/ntuple.jl)
function in `Base` uses this strategy to implement loop unrolling.

For benchmarks that compare these two implementations, see
[Manual vs. Recursive Unrolling](@ref).

## Interface API

The functions exported by this package can be used with any statically sized
iterators, as long as those iterators make appropriate use of the following
interface:

```@docs
generic_getindex
output_type_for_promotion
AmbiguousOutputType
NoOutputType
ConditionalOutputType
output_promote_rule
constructor_from_tuple
empty_output
StaticSequence
```

## How to Use the Interface

To unroll over a statically sized iterator of some user-defined type `T`, follow
these steps:
- Add a method for `getindex(::T, n)`, or for `generic_getindex(::T, n)` if
  `getindex` should not be defined for iterators of type `T`
- If every unrolled function that needs to construct an iterator when given an
  iterator of type `T` can return a `Tuple` instead, stop here
- Otherwise, to return a non-`Tuple` iterator whenever it is efficient to do so,
  follow these steps:
    - Add a method for `output_type_for_promotion(::T) = O`, where `O` can be
      `T`, a supertype of `T`, some other `Type`, or an `AmbiguousOutputType`
    - If an iterator whose output type is `O` can be used together with an
      iterator whose output type is `O′`, add a method for
      `output_promote_rule(O, O′)`
    - If `O` is a `NoOutputType`, stop here
    - Otherwise, to handle the unambiguous output type `U` that underlies `O`
      (where `U` is equivalent to `O` unless `O` is a `ConditionalOutputType`),
      follow these steps:
        - If an iterator of type `U` can be efficiently constructed from a
          `Tuple`, add a method for `constructor_from_tuple(U)`
        - Otherwise, for each of the following functions, add a method if it can
          be implemented to construct an iterator of type `U` without first
          storing the iterator's contents in a `Tuple`:
            - `empty_output(U)`
            - `unrolled_push_into(U, itr, item)`
            - `unrolled_append_into(U, itr1, itr2)`
            - `unrolled_take_into(U, itr, val_N)`
            - `unrolled_drop_into(U, itr, val_N)`
            - `unrolled_map_into(U, f, itr)`
            - `unrolled_accumulate_into(U, op, itr, init, transform)`

!!! note "Note"
    When a relevant method for the interface is not defined, unrolled functions
    will typically fall back to using `Tuple`s instead of other iterator types.
