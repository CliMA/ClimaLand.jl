# Jieko

[![CI](https://github.com/Roger-luo/Jieko.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Roger-luo/Jieko.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/Roger-luo/Jieko.jl/graph/badge.svg?token=8EIbN4OPo2)](https://codecov.io/gh/Roger-luo/Jieko.jl)
[![][docs-stable-img]][docs-stable-url]
[![][docs-dev-img]][docs-dev-url]

Documentation as interfaces. 接口 (Jiēkǒu) is the Chinese word for interfaces and APIs. This one works with the `public` keyword.

Julia uses docstrings to define interfaces. This is a flexible way of creating interfaces in a dynamic language, but also creates trouble for automation and tooling. Jieko is a package that provides a infrastructure for defining interfaces that works with DocStringExtension with precisely the signature of the interface.

## Installation

<p>
Jieko is a &nbsp;
    <a href="https://julialang.org">
        <img src="https://raw.githubusercontent.com/JuliaLang/julia-logo-graphics/master/images/julia.ico" width="16em">
        Julia Language
    </a>
    &nbsp; package. To install Jieko,
    please <a href="https://docs.julialang.org/en/v1/manual/getting-started/">open
    Julia's interactive session (known as REPL)</a> and press <kbd>]</kbd>
    key in the REPL to use the package mode, then type the following command
</p>

```julia
pkg> add Jieko
```

## Example

You only need to use the `@pub` macro and if you have [DocStringExtensions](https://github.com/JuliaDocs/DocStringExtensions.jl) setup, you can use the `DEF` stub to generate the interface definition in the docstring similar to the `SIGNATURES` or `TYPEDSIGNATURE` for methods.

```julia
using Jieko: @pub, DEF

"""
$DEF

my lovely interface
"""
@pub jieko(x::Real) = x
```

## Why Jieko?

Julia interfaces and public APIs are defined by documentation. There is no strict requirement
on which interface an object should implement, but rather the interface is defined by the
documentation. This is a flexible way of defining interfaces in a dynamic language, but it
also creates trouble for automation and tooling.

Existing approaches like [Interfaces](https://juliahub.com/ui/Packages/General/Interfaces) approach
the problem by defining the interface explicitly as an object. This is a good approach for people
want to define interfaces explicitly and in a more strict way. However, it becomes very verbose and not
working well with Julia's documenation system.

On the other hand, existing solution in [DocStringExtensions](https://juliahub.com/ui/Packages/General/DocStringExtensions)
fails to provide the precise interface definition in the docstring. For example, `SIGNATURES` will ignore
type annotations and only show the method name and argument names.

```julia
using DocStringExtensions: SIGNATURES

"""
$SIGNATURES

my lovely method
"""
doc_string_ext(x::Real) = x
```

they result in the following

```julia
help?> doc_string_ext
search: doc_string_ext

  doc_string_ext(x)
  

  my lovely method
```

On the other hand, `TYPEDSIGNATURE` can be too verbose, missing the type alias or messing up the return type.

```julia
using DocStringExtensions: TYPEDSIGNATURES

const MyAliasName = Int

"""
$TYPEDSIGNATURES

my lovely method
"""
doc_string_ext(x::Real)::Int = error("not implemented")

"""
$TYPEDSIGNATURES

my lovely method
"""
doc_string_ext(x::MyAliasName)::Complex = error("not implemented")
```

this results in the following

```julia
help?> doc_string_ext
search: doc_string_ext

  doc_string_ext(x::Real)
  

  my lovely method

  ─────────────────

  doc_string_ext(x::Int64)
  

  my lovely method
```

Using `@pub` and `DEF` fixes the problem as they actually record the precise interface definition.

```julia
using Jieko: @pub, DEF

const MyAliasName = Int

"""
$DEF
"""
@pub jieko(x::Real)::Int = error("not implemented")

"""
$DEF
"""
@pub jieko(x::MyAliasName)::Complex = error("not implemented")
```

this results in the following

```julia
help?> jieko
search: jieko Jieko

  jieko(x::Real) -> Int

  ─────────────────

  jieko(x::MyAliasName) -> Complex
```

In summary, the `@pub` macro from Jieko records the precise interface signature of your definition in the docstring, which can be used by tools to generate documentation or check the interface implementation.

See the [documentation](https://Roger-luo.github.io/Jieko.jl/dev/) for more details.

## License

MIT License

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://Roger-luo.github.io/Jieko.jl/dev/
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://Roger-luo.github.io/Jieko.jl/stable
