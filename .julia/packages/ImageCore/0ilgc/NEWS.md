# 0.10.0

## Breaking changes

- `clamp01` has been restricted to `AbstractGray` and `AbstractRGB`, as it did not make sense for colors like `HSV` (#181)
- Step 3 of [ColorVectorSpace's `abs2` transition has been completed](https://github.com/JuliaGraphics/ColorVectorSpace.jl#abs-and-abs2); users and developers should now replace any `ColorVectorSpace.Future.abs2` calls with `abs2` to prevent breakage when the final step completes.

## Other changes

- AbstractFFTs.jl is no longer a dependency; error hints are now used for `fft` operations on deliberately-unsupported `Colorant` types
- Graphics.jl is no longer a dependency
- PrecompileTools.jl is now used for precompilation

# 0.7.0

## Breaking changes

- Switch to Julia 0.7
- Because of changes to Julia's own `reinterpret`, ImageCore now
  defines and exports `reinterpretc` for reinterpreting
        numeric arrays <----> colorant arrays (#52)
- The ColorView and ChannelView types are deleted; their functionality
  is now handled via `reinterpret` (#52)

# 0.6.0

- Last minor version with Julia 0.6 support
