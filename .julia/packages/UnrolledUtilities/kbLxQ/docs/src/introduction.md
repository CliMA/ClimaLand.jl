```@setup inference_test
using UnrolledUtilities, InteractiveUtils, Test
```

```@setup fake_inference_test
macro code_warntype(_...) nothing end
macro code_llvm(_...) nothing end
```

```@raw html
<style>
summary {
  background-color: #3c5dcd;
  border-radius: 5px;
  color: white;
  cursor: pointer;
  list-style-position: outside;
  list-style-type: "⬇";
}
summary::after {
  content: "Click to show long output";
  margin-left: 15px;
}
details[open] summary {
  list-style-type: "⬆";
}
details[open] summary::after {
  content: none;
}
</style>
```

## Motivation for Loop Unrolling

Although the iteration utilities in `Base` and `Base.Iterators` are sufficiently
performant for most common use cases, those who choose to dive into the world of
low-level optimization will often discover
[type instabilities](https://docs.julialang.org/en/v1/manual/faq/#man-type-stability)
in unexpected situations. Here is a particularly simple example:

```@repl inference_test
Test.@inferred map(one, Tuple(1:31));
Test.@inferred map(one, Tuple(1:32));
```

This type instability is present in all `map`s over iterators with lengths
greater than 31, regardless of whether they are statically sized. As with most
type instabilities in Julia, this leads to memory allocations every time `map`
is called with sufficiently long iterators.

[`Test.@inferred`](https://docs.julialang.org/en/v1/stdlib/Test/#Test.@inferred)
is helpful for checking whether the return type of a function call is stable,
but looking directly at the generated [LLVM](https://llvm.org/docs/LangRef.html)
code reveals just how different the two function calls above are:

```@repl inference_test
@code_llvm debuginfo=:none map(one, Tuple(1:31))
```
```@raw html
<details><summary>
```
```@repl fake_inference_test
@code_llvm debuginfo=:none map(one, Tuple(1:32))
```
```@raw html
</summary>
```
```@repl inference_test
@code_llvm debuginfo=:none map(one, Tuple(1:32)) # hide
```
```@raw html
</details><br>
```

The type instability (and all of the resulting LLVM code complexity) in the
second function call can be eliminated by replacing `map` with `unrolled_map`:

```@repl inference_test
Test.@inferred unrolled_map(one, Tuple(1:32));
@code_llvm debuginfo=:none unrolled_map(one, Tuple(1:32))
```

The minimum iterator length for type instability is not always 32; for instance,
it can also be 14:

```@repl inference_test
first_11(itr) = itr[1:11]
Test.@inferred first_11(Tuple(1:13));
Test.@inferred first_11(Tuple(1:14));
```

!!! note "Note"
    ##### *Why is the function definition needed in this example?*

    On the first line of the example above, `[1:11]` is enclosed in a function
    so that it does not get evaluated in global scope. This turns the range
    `1:11` into a `Core.Const`, which the compiler can propagate into the call
    to `getindex` in order to infer the length of the result:

    ```@setup first_11_code_warntype
    using InteractiveUtils
    first_11(itr) = itr[1:11]
    ```

    ```@repl first_11_code_warntype
    @code_warntype first_11(Tuple(1:13))
    ```

    In contrast, running `Test.@inferred Tuple(1:13)[1:11]` would amount to
    checking whether the compiler can compute the result type of `getindex`
    given only the argument types `NTuple{13, Int64}` and `UnitRange{Int64}`,
    which it cannot do:

    ```@raw html
    <details><summary>
    ```
    ```@repl fake_inference_test
    @code_warntype Tuple(1:13)[1:11]
    ```
    ```@raw html
    </summary>
    ```
    ```@repl inference_test
    @code_warntype Tuple(1:13)[1:11] # hide
    ```
    ```@raw html
    </details><br>
    ```

Although `itr[1:10]` is always inferrable when `itr` is a `Tuple`, `itr[1:11]`
has a type instability whenever `itr` contains more than 13 items. More
generally, `itr[1:N]` seems to be unstable for all `N > 10` whenever `itr`
contains more than `N + 2` items. This type instability can be fixed by
replacing `getindex` with `unrolled_take`:

```@repl inference_test
unrolled_first_11(itr) = unrolled_take(itr, Val(11))
Test.@inferred unrolled_first_11(Tuple(1:14));
```

Even when the final result of a function is inferred, there can be intermediate
steps in the function with type instabilities that trigger allocations:

```@repl inference_test
function add_lengths(itr)
    length_sum = 0
    for n in 1:length(itr)
        length_sum += length(itr[n])
    end
end
Test.@inferred add_lengths(((1, 2), (1, 2, 3)))
@allocated add_lengths(((1, 2), (1, 2, 3)))
@code_warntype add_lengths(((1, 2), (1, 2, 3)))
```

The output of `@code_warntype` is quite cluttered, but the most important detail
here is that the call to `getindex` does not get inferred because it can result
in either a `Tuple` of length 2 or a `Tuple` of length 3. This type instability
can be fixed by replacing `getindex` with `unrolled_applyat`:

```@repl inference_test
function unrolled_add_lengths(itr)
    length_sum = 0
    for n in 1:length(itr)
        length_sum += unrolled_applyat(length, n, itr)
    end
end
unrolled_add_lengths(((1, 2), (1, 2, 3))) # hide
@allocated unrolled_add_lengths(((1, 2), (1, 2, 3)))
@code_warntype unrolled_add_lengths(((1, 2), (1, 2, 3)))
```

For a detailed breakdown of when the tools provided by this package can improve
performance, see the [User Guide](user_guide.md).

## What Does Loop Unrolling Do

When a loop over `N` indices is unrolled, it gets compiled into `N` lines of
LLVM code, where each line has a constant (`Core.Const`) index. For example, an
unrolled loop that prints every integer from 1 to 33 is compiled into the
following:

```@raw html
<details><summary>
```
```@repl fake_inference_test
@code_llvm debuginfo=:none unrolled_foreach(println, Tuple(1:33))
```
```@raw html
</summary>
```
```@repl inference_test
@code_llvm debuginfo=:none unrolled_foreach(println, Tuple(1:33)) # hide
```
```@raw html
</details><br>
```

This LLVM code consists of 33 `getelementptr` instructions (each of which
extracts a value from a `Tuple` at a particular index), 33 `load` instructions,
and 33 `call` instructions (each of which switches execution to `println`).
Every `getelementptr` instruction has a constant index between 0 and 32; in more
complex examples where the `call` instructions get inlined, this constant index
can be propagated into the LLVM code of the function being called. On the other
hand, here is the LLVM code for the non-unrolled version of this loop:

```@repl inference_test
@code_llvm debuginfo=:none foreach(println, Tuple(1:33))
```

This LLVM code has a `load` instruction to get the first value and a
`getelementptr` instruction with a non-constant integer index to get all other
values. It also has conditional jump instructions for checking whether the last
index has been reached after each `load` and `getelementptr` instruction.

## Downsides of Loop Unrolling

Given the performance benefits of loop unrolling, it might seem at first that
the standard library needs more of it. However, the standard library is not just
meant for writing high-performance code with statically sized iterators—many of
its use cases involve code that is only executed once or several times. In such
cases, most of the execution time is required for compilation, and minimizing
run time makes no practical difference. Although unrolled functions can
occasionally be faster to compile than non-unrolled functions, they are
typically slower to compile, which means that using them instead of standard
library functions can often increase total execution time:

```@repl inference_test
tup32 = ntuple(Returns((1, 2)), 32);
@elapsed map(first, tup32)
@elapsed unrolled_map(first, tup32)
```

The increase in compilation time is usually no more than a factor of 5 for small
iterators, but it grows as iterator length increases:

```@repl inference_test
tup320 = ntuple(Returns((1, 2)), 320);
@elapsed map(first, tup320)
@elapsed unrolled_map(first, tup320)
```

Loop unrolling can also theoretically increase the run time of a function in
addition to its compilation time, since unrolled assembly code requires more
space and takes longer to load than non-unrolled code. In practice, though, the
constant propagation enabled by unrolling usually compensates for this slowdown.

So, when type instabilities and memory allocations need to be removed
([as is required for static compilation](https://github.com/brenhinkeller/StaticTools.jl#limitations))
and the cost to total execution time is more or less irrelevant, using unrolled
functions is probably worthwhile. Otherwise, if a significant increase in
compilation time (and potentially also run time) needs to be avoided, using
standard library functions might be a better option.

It is usually a good idea to compare the performance of unrolled code against
non-unrolled code before settling on a particular design. Many examples of such
comparisons can be found in the [tables of benchmarks](comparison_tables.md)
that are automatically generated for this package.
