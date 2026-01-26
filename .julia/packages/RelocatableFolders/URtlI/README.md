# RelocatableFolders.jl

An alternative to the `@__DIR__` macro. Packages that wish to reference paths
in their project directory run into issues with relocatability when used in
conjunction with `PackageCompiler`. The `@path` macro provided by this package
overcomes this limitation. See [here][pkgcompiler] and [here][julia-issue] for
further details.

[pkgcompiler]: https://julialang.github.io/PackageCompiler.jl/stable/apps.html#relocatability
[julia-issue]: https://github.com/JuliaLang/julia/issues/38696

## Usage

The package provides one export, the `@path` macro. It can be used to replace
`@__DIR__` in the following way:

```julia
module MyPackage

using RelocatableFolders

# const ASSETS = joinpath(@__DIR__, "../assets")
const ASSETS = @path joinpath(@__DIR__, "../assets")

end
```

At *runtime* the path stored in `ASSETS` will get resolved to either the
original path, if it still exists, or to an automatically generated
scratchspace containing the same folder and file structure as the original.

## Limitations

This macro should only be used for reasonably small file or folder sizes. If
there are very large files then it is better to make use of Julia's `Artifact`
system instead.

Building new paths from, for example, `ASSETS` in the above example will return
a `String` containing the resolved path rather than a `Path` object. Doing this
at the module-level will result in hardcoded paths that will run into
relocatability issues as discussed above. Always create a new `@path` for
each resource you wish to reference rather than building them in parts, e.g.

```julia
module MyPackage

using RelocatableFolders

const ASSETS = @path joinpath(@__DIR__, "../assets")
const SUBDIR = @path joinpath(ASSETS, "subdir")
const FILE = @path joinpath(ASSETS, "file.txt")

end
```

## Internals

At compile-time the `@path` macro will read in all the files contained in the
referenced path and store them and their paths. The returned object is a
`Path <: AbstractString`. Whenever a `Path` is passed to a function
expecting an `AbstractString` (such as `readdir`) it will be converted to a
`String` by looking up the stored path and returning that. When no path exists
(the source tree no longer exists) then the contents of the files that were
read at compile-time are written to a `Scratch` scratchspace and that path is
returned instead.

## Alternatives

The alternative approach is to use the `Artifacts` system to distribute the
required files, which is a reasonably heavyweight solution for a simple
collection of source-controlled files. It is recommended that users look to use
artifacts when that file sizes are reasonably large and only use this package
when the distributed files are small.
