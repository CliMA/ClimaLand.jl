```@setup inference_test
using UnrolledUtilities
```

#  UnrolledUtilities.jl

A toolkit for low-level optimization of Julia code in which iterator sizes are
known during compilation.

This package can be used with all *statically sized* iterators (`Tuple`s,
`NamedTuple`s, [`StaticArray`s](https://github.com/JuliaArrays/StaticArrays.jl),
etc.), including ones that are very long or ones that have elements of different
types, both of which are cases that Julia's standard library often handles
inefficiently. For example, the standard libary function `in` performs worse
than this package's `unrolled_in` for `Tuple`s with elements of different types:

```@repl inference_test
nonuniform_itr = ((1, 2), (1, 2, 3));
() in nonuniform_itr # hide
@allocated () in nonuniform_itr
unrolled_in((), nonuniform_itr) # hide
@allocated unrolled_in((), nonuniform_itr)
```

The [loop unrolling](https://en.wikipedia.org/wiki/Loop_unrolling) automatically
performed by this package offers the following benefits for statically sized
iterators:
- better support for *static compilation*
  - compilation of [executables](https://github.com/tshort/StaticCompiler.jl)
  - compilation of [GPU kernels](https://github.com/JuliaGPU/CUDA.jl)
- better performance (usually)
  - reduced run times
  - reduced memory footprints while code is running
- better compilation efficiency (occasionally)
  - reduced compilation times
  - reduced memory footprints while code is compiling

To find out more about loop unrolling and when it is useful, see the
[Introduction](introduction.md).

## Package Features

This package exports a number of analogues to functions from `Base` and
`Base.Iterators`, each of which has been optimized for statically sized
iterators (in terms of both performance and compilation time):
- `unrolled_push(itr, item)`—similar to `push!`, but non-mutating
- `unrolled_append(itr, itrs...)`—similar to `append!`, but non-mutating
- `unrolled_prepend(itr, itrs...)`—similar to `prepend!`, but non-mutating
- `unrolled_map(f, itrs...)`—similar to `map`
- `unrolled_any([f], itr)`—similar to `any`
- `unrolled_all([f], itr)`—similar to `all`
- `unrolled_foreach(f, itrs...)`—similar to `foreach`
- `unrolled_reduce(op, itr; [init])`—similar to `reduce` (i.e., `foldl`)
- `unrolled_mapreduce(f, op, itrs...; [init])`—similar to `mapreduce` (i.e.,
  `mapfoldl`)
- `unrolled_accumulate(op, itr; [init])`—similar to `accumulate`
- `unrolled_in(item, itr)`—similar to `in`
- `unrolled_unique([f], itr)`—similar to `unique`
- `unrolled_allunique([f], itr)`—similar to `allunique`
- `unrolled_allequal([f], itr)`—similar to `allequal`
- `unrolled_sum([f], itr; [init])`—similar to `sum`, but with `init = 0` when
  `itr` is empty
- `unrolled_prod([f], itr; [init])`—similar to `prod`, but with `init = 1` when
  `itr` is empty
- `unrolled_cumsum([f], itr)`—similar to `cumsum`, but with an optional `f`
- `unrolled_cumprod([f], itr)`—similar to `cumprod`, but with an optional `f`
- `unrolled_count([f], itr)`—similar to `count`
- `unrolled_maximum([f], itr)`—similar to `maximum`
- `unrolled_minimum([f], itr)`—similar to `minimum`
- `unrolled_extrema([f], itr)`—similar to `extrema`
- `unrolled_findmax([f], itr)`—similar to `findmax`
- `unrolled_findmin([f], itr)`—similar to `findmin`
- `unrolled_argmax([f], itr)`—similar to `argmax`
- `unrolled_argmin([f], itr)`—similar to `argmin`
- `unrolled_findfirst([f], itr)`—similar to `findfirst`
- `unrolled_findlast([f], itr)`—similar to `findlast`
- `unrolled_filter(f, itr)`—similar to `filter`
- `unrolled_flatten(itr)`—similar to `Iterators.flatten`
- `unrolled_flatmap(f, itrs...)`—similar to `Iterators.flatmap`
- `unrolled_product(itrs...)`—similar to `Iterators.product`
- `unrolled_cycle(itr, ::Val{N})`—similar to `Iterators.cycle`, but with a
  static value of `N`
- `unrolled_partition(itr, ::Val{N})`—similar to `Iterators.partition`, but with
  a static value of `N`
- `unrolled_take(itr, ::Val{N})`—similar to `Iterators.take` (i.e., `itr[1:N]`),
  but with a static value of `N`
- `unrolled_drop(itr, ::Val{N})`—similar to `Iterators.drop` (i.e.,
  `itr[(N + 1):end]`), but with a static value of `N`

In addition, this package exports several functions that do not have analogues
in `Base` or `Base.Iterators`:
- `unrolled_applyat(f, n, itrs...)`—similar to `f(itrs[1][n], itrs[2][n], ...)`
- `unrolled_argfirst(f, itr)`—similar to `itr[findfirst(f, itr)]`
- `unrolled_arglast(f, itr)`—similar to `itr[findlast(f, itr)]`
- `unrolled_split(f, itr)`—similar to `(filter(f, itr), filter(!f, itr))`, but
  without duplicate calls to `f`

These unrolled functions are compatible with the following types of iterators:
- statically sized iterators from `Base` (e.g., `Tuple` and `NamedTuple`)
- statically sized iterators from `StaticArrays` (e.g., `SVector` and `MVector`)
- lazy iterators from `Base` (e.g., the results of generator expressions,
  `Iterators.map`, `Iterators.reverse`, `enumerate`, and `zip`) that are used as
  wrappers for statically sized iterators

They are also compatible with two new types of statically sized iterators
exported by this package:
- `StaticOneTo`—similar to `Base.OneTo`
- `StaticBitVector`—similar to `BitVector`

See the [User Guide](@ref "When to Use StaticOneTo and StaticBitVector") for
additional information about these new types of iterators.

See the [Developer Guide](@ref "How to Use the Interface") to learn how
user-defined iterator types can be made compatible with unrolled functions.
