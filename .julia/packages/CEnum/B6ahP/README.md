# CEnum

[![Build Status](https://github.com/JuliaInterop/CEnum.jl/workflows/CI/badge.svg)](https://github.com/JuliaInterop/CEnum.jl/actions)
[![Codecov](https://codecov.io/gh/JuliaInterop/CEnum.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaInterop/CEnum.jl)
[![pkgeval](https://juliahub.com/docs/CEnum/pkgeval.svg)](https://juliahub.com/ui/Packages/CEnum/0gyUJ)
[![version](https://juliahub.com/docs/CEnum/version.svg)](https://juliahub.com/ui/Packages/CEnum/0gyUJ)
[![deps](https://juliahub.com/docs/CEnum/deps.svg)](https://juliahub.com/ui/Packages/CEnum/0gyUJ?t=2)
[![Genie Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/CEnum)](https://pkgs.genieframework.com?packages=CEnum)

This package provides a C-compatible enum for Julia.

```julia
julia> @enum Foo a = 1 b = 2 c = 1
ERROR: LoadError: ArgumentError: values for Enum Foo are not unique
Stacktrace:
 [1] @enum(::LineNumberNode, ::Module, ::Any, ::Vararg{Any,N} where N) at ./Enums.jl:128
in expression starting at REPL[12]:1

julia> using CEnum

julia> @cenum(Bar, d = 1, e = 2, f = 1)

julia> d == f
true
```

## Credit
The original version is written by [SimonDanisch](https://github.com/SimonDanisch) in [this PR](https://github.com/JuliaInterop/Clang.jl/pull/162).
