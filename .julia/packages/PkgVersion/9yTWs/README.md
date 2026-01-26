# PkgVersion.jl

[![Build Status][gha-img]][gha-url]    [![Coverage Status][codecov-img]][codecov-url]

Provide macros to access fields `version`, `uuid`, `authors` in `Project.toml at compile time.

## Usage

```julia
module MyModule
using Tar
using PkgVersion

const VERSION = PkgVersion.@Version 0
const UUID = PkgVersion.@Uuid 
const AUTHOR = PkgVersion.@Author "unknown@nowhere"

const VERSION_TAR = PkgVersion.Version(Tar)  

end
```

## Notes

File `Project.toml` must be readable in package directory.

`Author` returns the first string of field `authors`.

[gha-img]: https://github.com/KlausC/PkgVersion.jl/workflows/CI/badge.svg
[gha-url]: https://github.com/KlausC/PkgVersion.jl/actions?query=workflow%3ACI

[codecov-img]: https://codecov.io/gh/KlausC/PkgVersion.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/KlausC/PkgVersion.jl
