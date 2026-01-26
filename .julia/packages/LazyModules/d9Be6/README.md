# LazyModules

> No, no, not now

This package provides package developers an alternative option to delay package loading until used.
If some dependency is not used, then users don't need to pay for its latency.

This package is not panacea, it only works for a very limited set of use cases. This package is only
for (maybe experienced) package authors. End-users is not recommended to use this package directly.

Be careful about the world age issue!

## Examples

The eager mode is checked by setting environment variable `LazyModules_lazyload` as `"false"`.
All the functionalities will still be the same except that lazy mode delays a significant portion
of packages loading to their first usage -- if ever get called.

| Package                 | latency (eager)    | latency (lazy) | Julia   |
| ------------------      | -----------------  | -------------- | ------  |
| [ImageIO.jl] v0.6.3     |  1.049s            |  0.115s        | v1.7.2  |

The latency is measured by `@time using XXX` in a new Julia REPL after full precompilation.  Each
row is measured on the same machine.

## Syntax

- `@lazy import Foo = "<UUID>"` ✅
- `@lazy import Foo as LazyFoo = "<UUID>"` ✅ (Julia 1.6+)
- `@lazy import Foo` ❌
- `@lazy using Foo` ❌

## The lazy Plots story

Assume that you've built a fantastic package `examples/MyPkg` with some built-in plot functions:

```julia
module MyPkg

export generate_data, draw_figure
import Plots

generate_data(n) = sin.(range(start=0, stop=5, length=n) .+ 0.1.*rand(n))
draw_figure(data) = plot(data, title="MyPkg Plot")

end
```

Normally, you spend quite a long time on loading the package because `Plots` is heavy:

```julia
(@v1.7) pkg> activate examples/MyPkg
  Activating project at `~/Documents/Julia/LazyModules.jl/examples/MyPkg`

(MyPkg) pkg> instantiate
Precompiling project...
  1 dependency successfully precompiled in 36 seconds (133 already precompiled)

julia> @time using MyPkg # 💤
  2.857596 seconds (9.81 M allocations: 670.470 MiB, 8.53% gc time, 19.95% compilation time)

julia> x = @time generate_data(100); # 🚀
  0.000006 seconds (2 allocations: 1.750 KiB)

julia> @time draw_figure(x) # 💤
1.608146 seconds (4.00 M allocations: 223.266 MiB, 2.83% gc time, 99.74% compilation time)
```

If Plots is the needed feature to `MyPkg`, then the latency is what I need to pay for, which is
okay. **BUT**, from time to time, I might just generate the data and save it to disk, **without
plotting the figure at all!** Then why should I still wait for the `Plots` loading?

This is where `LazyModules` can become useful: it delays the loading of heavy packages such as
`Plots` to its first call. By doing this, we don't need to wait for it if we don't use the `Plots`
functionalities.

What you need to do, is to change the package code a bit (`examples/MyLazyPkg`):

```diff
module MyLazyPkg

export generate_data, draw_figure
+using LazyModules
-import Plots
+@lazy import Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80"

generate_data(n) = sin.(range(start=0, stop=5, length=n) .+ 0.1.*rand(n))
draw_figure(data) = Plots.plot(data, title="MyPkg Plot")

end
```

By doing this, if the users don't use `draw_figure` feature, then they don't need to load `Plots` at
all, which makes package loading significantly faster:

```julia
(@v1.7) pkg> activate examples/MyLazyPkg
  Activating project at `~/Documents/Julia/LazyModules.jl/examples/MyPkg`

(MyLazyPkg) pkg> instantiate
Precompiling project...
  1 dependency successfully precompiled in 36 seconds (133 already precompiled)

julia> @time using MyLazyPkg # 🚀🚀🚀🚀🚀
  0.053273 seconds (154.16 k allocations: 8.423 MiB, 97.62% compilation time)

julia> x = @time generate_data(100); # 🚀
  0.000006 seconds (2 allocations: 1.750 KiB)
```

The actual loading of `Plots` is delayed to the first `draw_figure` call:

```julia
julia> @time draw_figure(x) # 💤💤
  4.454738 seconds (13.82 M allocations: 897.071 MiB, 8.81% gc time, 49.97% compilation time)
```

Here `4.4` seconds is approximately `2.8` (Plots loading time) plus `1.6` (time to first plot). For
this reason, if a functionality is really necessary and widely used by almost everyone, then this
LazyModules package won't be helpful at all.

## What is a LazyModule

`LazyModule` is not a `Module`; it is indeed, a struct that overrides `getproperty`.

```julia
julia> using LazyModules

julia> @lazy import SparseArrays="2f01184e-e22b-5df5-ae63-d93ebab69eaf"
LazyModule(SparseArrays)

julia> SparseArrays.sprand(10, 10, 0.3) # triggers the loading
10×10 SparseArrays.SparseMatrixCSC{Float64, Int64} with 40 stored entries:
...
```

Package is loaded whenever there's a `getproperty` call, e.g., `SparseArrays.sprand` as shown above.

## World-age issue

The simplest example to trigger the world age issue is perhaps the following:

```julia
julia> using LazyModules

julia> @lazy import ImageCore = "a09fc81d-aa75-5fe9-8630-4744c3626534"
LazyModule(ImageCore)

julia> function foo()
           c = ImageCore.RGB(0.0, 0.0, 0.0)
           return c .* 3
       end
foo (generic function with 1 method)

julia> foo()
ERROR: MethodError: no method matching length(::ColorTypes.RGB{Float64})
The applicable method may be too new: running in world age 31343, while current world is 31370.
...

julia> foo()
RGB{Float64}(0.0,0.0,0.0)
```

Here we can see that:

- at first `foo()` call, it triggers the world-age issue
- at the second call, it is working okay

This happens because when you first call `foo()`, the `length` method required by `*` is not yet
available (to the current world age). When the `ImageCore.RGB` triggers the package loading of
`ImageCore`, which again triggers the recompilation of many methods (in a new world age). But still,
`*` from the old world age can't see the `length` method in the new world age. Things changed at the
second call, where `foo()` gets recompiled in the new world age.

There are commonly two ways to work around the world-age issue:

The first workaround is to use `invokelatest` whenever world-age issue occurs.
But this has some overhead due to the dynamic dispatch.

```julia
julia> using LazyModules

julia> @lazy import ImageCore
LazyModule(ImageCore)

julia> function foo()
           c = ImageCore.RGB(0.0, 0.0, 0.0)
           return Base.invokelatest(*, c, 3)
       end
foo (generic function with 1 method)

julia> foo()
RGB{Float64}(0.0,0.0,0.0)
```

The second workaround is to load the "core" packages eagerly so that we don't need to process
"alien" types. For instance, `RGB` and its arithmetic are provided by `Colors` and
`ColorVectorSpace`:

```julia
julia> using Colors, ColorVectorSpace

julia> using LazyModules

julia> @lazy import ImageCore = "a09fc81d-aa75-5fe9-8630-4744c3626534"
LazyModule(ImageCore)

julia> function foo()
           c = ImageCore.RGB(0.0, 0.0, 0.0)
           return c * 3
       end
foo (generic function with 1 method)

julia> foo()
RGB{Float64}(0.0,0.0,0.0)
```

The world-age issue is exactly the reason why this package should not be used by users directly.

As you might already notice, the main cause of world age issue is because some packages are not
loaded eagerly. This package provides a helper function to require user codes eagerly load the
requested packages. For instance,

```julia
julia> @lazy import ImageCore as LazyImageCore = "a09fc81d-aa75-5fe9-8630-4744c3626534"
LazyModule(ImageCore)

julia> function foo()
           LazyModules.require(LazyImageCore)
           c = LazyImageCore.RGB(0.0, 0.0, 0.0)
           return c * 3
       end
foo (generic function with 1 method)

julia> foo()
ERROR: ImageCore is required to be loaded first, maybe `using ImageCore` or `import ImageCore` and try again.
...

julia> using ImageCore # now let's explicitly load it

julia> foo() # issue gone
RGB{Float64}(0.0,0.0,0.0)
```

By doing this, your users can get an informative error messages rather than arbitrary non-sense
errors due to world-age issue. This can be used to "lock" some features unless the users explicitly
load necessary dependencies.

## Overhead

The overhead is about ~50ns in Intel i9-12900K due to the dynamic dispatch via `invokelatest`. Thus
you should not use this package for very trivial functions.

```julia
julia> @lazy import ImageCore as LazyImageCore = "a09fc81d-aa75-5fe9-8630-4744c3626534"
LazyModule(ImageCore)

julia> import ImageCore

julia> @btime zero(LazyImageCore.RGB)
  49.004 ns (0 allocations: 0 bytes)
RGB{N0f8}(0.0,0.0,0.0)

julia> @btime zero(ImageCore.RGB)
  0.012 ns (0 allocations: 0 bytes)
RGB{N0f8}(0.0,0.0,0.0)
```

## Disable lazy mode

If you want to eagerly load every dependency, you can start your Julia process with
environment variable `LazyModules_lazyload` as `"false"`. That way `@lazy import` command
will degenerated to a normal `import` command. For instance, in bash shell, you could do

```console
jc@ubuntu:~$ LazyModules_lazyload=false julia
```

If set successfully, this will print an informative message:

```julia
[ Info: disable lazy package loading: environment variable `LazyModules_lazyload=false` is detected
```

[ImageIO.jl]: https://github.com/JuliaIO/ImageIO.jl
