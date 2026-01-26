```@meta
CurrentModule = UnrolledUtilities
```

```@setup inference_test
using UnrolledUtilities, InteractiveUtils, Test
nonuniform_itr = ((1, 2), (1, 2, 3));
nested_itr_of_depth_2 = ((1, 2), (3, 4));
nested_itr_of_depth_3 = (((1,), (2,)), ((3,), (4,)));
```

# When to Use UnrolledUtilities

The functions and types exported by this package tend to perform better than
their counterparts from `Base` and `Base.Iterators` in the scenarios listed
below. Additional examples and more precise measurements can be found in the
automatically generated [tables of benchmarks](comparison_tables.md).

##### Outline:

```@contents
Pages = ["user_guide.md"]
Depth = 2:3
```

## When to Use Unrolled Functions

### Long iterators

- `map` has an unstable return type for iterators with lengths greater than 32:

  ```@repl inference_test
  Test.@inferred map(one, Tuple(1:31));
  Test.@inferred map(one, Tuple(1:32));
  Test.@inferred unrolled_map(one, Tuple(1:32));
  ```

- `getindex` has an unstable return type for `Core.Const` slices of length
  `N > 10` from iterators with lengths greater than `N + 2`:

  ```@repl inference_test
  first_11(itr) = itr[1:11]
  Test.@inferred first_11(Tuple(1:13));
  Test.@inferred first_11(Tuple(1:14));
  unrolled_first_11(itr) = unrolled_take(itr, Val(11))
  Test.@inferred unrolled_first_11(Tuple(1:14));
  ```

- For benchmarks that indicate performance improvements when using unrolled
  functions with long iterators, see [Isolated Unrolled Functions](@ref)

### Iterators with elements of different types

- `in`, `any`, `all`, `foreach`, and other functions in `Base` have intermediate
  type instabilities that trigger allocations for nonuniform iterators:

  ```@repl inference_test
  nonuniform_itr = ((1, 2), (1, 2, 3));
  () in nonuniform_itr # hide
  @allocated () in nonuniform_itr
  unrolled_in((), nonuniform_itr) # hide
  @allocated unrolled_in((), nonuniform_itr)
  any(isempty, nonuniform_itr) # hide
  @allocated any(isempty, nonuniform_itr)
  unrolled_any(isempty, nonuniform_itr) # hide
  @allocated unrolled_any(isempty, nonuniform_itr)
  ```

- `getindex` has an unstable return type for nonuniform iterators when given
  non-constant (i.e., not `Core.Const`) indices, which can lead to intermediate
  type instabilities that trigger allocations:

  ```@repl inference_test
  function add_lengths(itr)
      length_sum = 0
      for n in 1:length(itr)
          length_sum += length(itr[n])
      end
  end
  add_lengths(nonuniform_itr) # hide
  @allocated add_lengths(nonuniform_itr)
  function unrolled_add_lengths(itr)
      length_sum = 0
      for n in 1:length(itr)
          length_sum += unrolled_applyat(length, n, itr)
      end
  end
  unrolled_add_lengths(nonuniform_itr) # hide
  @allocated unrolled_add_lengths(nonuniform_itr)
  ```

  !!! note "Note"
      ##### *How can `unrolled_applyat` be stable if `n` isn't a `Core.Const`?*

      For the example of `add_lengths`, the compiler must infer the return
      type of `itr[::Int64]` before it can compile the call to `length`.
      Since this return type depends on the index `n`, the compiler needs to
      insert a runtime lookup into the method table that determines which
      method of `length` to call, `length(::Tuple{Int64, Int64})` or
      `length(::Tuple{Int64, Int64, Int64})`, and this triggers allocations.

      For the example of `unrolled_add_lengths`, the compiler instead infers
      the return types of `itr[::Core.Const(1)]`, `itr[::Core.Const(2)]`,
      and so on for every index into `itr`. Then, it compiles a call to
      `length` for each of these return types, and it inserts a runtime
      [switch instruction](https://llvm.org/docs/LangRef.html#switch-instruction)
      that determines which result of `length` to return for a particular
      value of `n`. As long as `length` itself only returns one type (in this
      case, `Int64`), this ensures that `unrolled_add_lengths` has no
      intermediate type instabilities.
      
      In other words, `unrolled_applyat` combines multiple methods for `length`
      and `getindex` into a single method, replacing the inefficient method
      table lookup that switches between them with a simpler switch instruction.

  !!! tip "Tip"
      ##### *When should `getindex` be replaced with `unrolled_applyat`?*

      The specific example above can be simplified by using `unrolled_sum`,
      instead of using a `for`-loop in conjunction with `unrolled_applyat`:

      ```@repl inference_test
      nonuniform_itr = ((1, 2), (1, 2, 3));
      unrolled_sum(length, nonuniform_itr) # hide
      @allocated unrolled_sum(length, nonuniform_itr)
      ```

      However, there are often situations in which loops cannot be replaced with
      function calls, such as when the loops are parallelized over multiple
      threads. For example, since CUDA is unable to compile kernels with type
      instabilities, something like the switch instruction in `unrolled_applyat`
      is *required* in order to parallelize over nonuniform iterators on GPUs.

- For benchmarks that indicate performance improvements when using unrolled
  functions with nonuniform iterators, see [Isolated Unrolled Functions](@ref)
  and [Nested Unrolled Functions](@ref)

### Reductions with intermediate values of different types

- `reduce` and `accumulate` have unstable return types when the return type of
  the reduction operator is not constant, but only for iterators with lengths
  greater than 32:

  ```@repl inference_test
  Test.@inferred reduce(tuple, Tuple(1:32));
  Test.@inferred reduce(tuple, Tuple(1:33));
  Test.@inferred unrolled_reduce(tuple, Tuple(1:33));
  ```

- For benchmarks that indicate performance improvements when using unrolled
  functions with nonuniform reductions, see [Isolated Unrolled Functions](@ref)

### Functions with recursion during compilation

- Any function that recursively calls itself with different types of arguments
  can trigger a compilation heuristic called the "recursion limit", which causes
  the compiler to generate code with type instabilities and allocations
  
  - Before Julia 1.11, the default recursion limit applied to every function
    with at least 3 levels of recursion during compilation, including relatively
    simple functions like an analogue of `length` for nested `Tuple`s:

    ```@repl inference_test
    nested_itr_of_depth_2 = ((1, 2), (3, 4));
    nested_itr_of_depth_3 = (((1,), (2,)), ((3,), (4,)));
    recursive_length(itr) =
        eltype(itr) <: Tuple ? sum(recursive_length, itr) : length(itr)
    Test.@inferred recursive_length(nested_itr_of_depth_2);
    Test.@inferred recursive_length(nested_itr_of_depth_3);
    ```
  
  - As of Julia 1.11, the default recursion limit is no longer triggered for
    this simple example, but still applies to more complex recursive functions:

    ```@repl inference_test
    recursive_sum_with_logs(itr) =
        eltype(itr) <: Tuple ?
        sum(recursive_sum_with_logs, itr) +
        sum(log ∘ recursive_sum_with_logs, itr) :
        sum(itr)
    Test.@inferred recursive_sum_with_logs(nested_itr_of_depth_2);
    Test.@inferred recursive_sum_with_logs(nested_itr_of_depth_3);
    unrolled_recursive_sum_with_logs(itr) =
        eltype(itr) <: Tuple ?
        unrolled_sum(unrolled_recursive_sum_with_logs, itr) +
        unrolled_sum(log ∘ unrolled_recursive_sum_with_logs, itr) :
        unrolled_sum(itr)
    Test.@inferred unrolled_recursive_sum_with_logs(nested_itr_of_depth_3);
    ```

  !!! note "Note"
      ##### *How can the default recursion limit be avoided?*

      The recursion limit of any function `f` can be disabled by evaluating the
      following code in the module where `f` is defined:

      ```julia
      @static if hasfield(Method, :recursion_relation)
          for method in methods(f)
              method.recursion_relation = Returns(true)
          end
      end
      ```

      However, the recursion limits of functions defined in `Base` cannot be
      modified from any module outside of `Base`. Since the default limit
      applies to all functions from `Base`, they can have unstable return types
      whenever they need more than 2 levels of recursion for compilation, even
      if the user-defined functions passed to them have no recursion limits. The
      only way to avoid the default limit of a function from `Base` is to
      replace it with a function whose limit has been disabled (such as an
      analogous function from `UnrolledUtilities`).

- For benchmarks that indicate performance improvements when using unrolled
  functions with recursive operations, see [Recursive Unrolled Functions](@ref)

## When to Use `StaticOneTo` and `StaticBitVector`

### Iterators of `Int`s from 1 to `N`

```@docs
StaticOneTo
```

If an iterator only contains the integers from 1 to `N ≥ 0`, it is possible to
provide the compiler with the values in the iterator in addition to their types
by using a `StaticOneTo`, as opposed to a `Tuple` or something similar. This
can allow the compiler to fully optimize out code that depends on those values,
essentially moving the code's execution from run time to compilation time:

```@repl inference_test
@code_llvm debuginfo=:none mapreduce(abs2, +, (1, 2, 3))
@code_llvm debuginfo=:none mapreduce(abs2, +, StaticOneTo(3))
```

Standard library functions can sometimes take advantage of this optimization,
but for most non-trivial operations it is necessary to use unrolled functions:

```@repl inference_test
@code_llvm debuginfo=:none mapreduce(log, +, StaticOneTo(3))
@code_llvm debuginfo=:none unrolled_mapreduce(log, +, StaticOneTo(3))
```

For benchmarks that indicate performance improvements when using `StaticOneTo`s,
see [Very Long Iterators](@ref).

!!! note "Note"
    ##### *Can the compiler infer iterator values in other scenarios?*

    The compiler can usually infer the values of iterators that only contain
    [singletons](https://docs.julialang.org/en/v1/manual/types/#man-singleton-types)
    when they are accessed using `Core.Const` indices, but this is not possible
    for non-singletons (e.g., integers) unless some special type of iterator is
    used (e.g., a `StaticOneTo`).

### Long iterators of `Bool`s that get modified across loop iterations

```@docs
StaticBitVector
```

Loops in Julia often allocate memory when a value larger than 32 bytes in size
is modified across loop iterations (regardless of whether the loops are unrolled
or not). Since `Bool`s are represented by bytes, this limits certain types of
loops to modifying [bitmasks](https://en.wikipedia.org/wiki/Mask_(computing)) of
no more than 32 `Bool`s in order to avoid allocations. Unlike an iterator of
`Bool`s, though, a `StaticBitVector` stores 8 bits in every byte, which makes it
possible to modify up to 256 bits at a time in loops without any allocations:

```@repl inference_test
random_bit_flips(itr) = reduce(
    (itr′, i) -> Base.setindex(itr′, !itr′[rand(1:i)], i),
    1:length(itr);
    init = itr,
)
@allocated random_bit_flips(ntuple(Returns(true), Val(32))) # hide
@allocated random_bit_flips(ntuple(Returns(true), Val(32)))
@allocated random_bit_flips(ntuple(Returns(true), Val(33))) # hide
@allocated random_bit_flips(ntuple(Returns(true), Val(33)))
@allocated random_bit_flips(StaticBitVector{256}(true)) # hide
@allocated random_bit_flips(StaticBitVector{256}(true))
```

As with `StaticOneTo`s, standard library functions can occasionally optimize
`StaticBitVector`s as well as unrolled functions, but most complex use cases
require unrolled functions.

For benchmarks that indicate performance improvements when using long
`StaticBitVector`s that get modified across loop iterations, see
[Nested Unrolled Closures](@ref).
