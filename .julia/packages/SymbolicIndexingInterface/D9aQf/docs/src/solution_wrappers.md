# Defining Solution Wrapper Fallbacks

The simplest case is when the type contains an object that already implements the interface.
All its methods can simply be forwarded to that object. To do so, SymbolicIndexingInterface.jl
provides the [`symbolic_container`](@ref) method. For example,

```julia
struct MySolutionWrapper{T<:SciMLBase.AbstractTimeseriesSolution}
  sol::T
  # other properties...
end

symbolic_container(sys::MySolutionWrapper) = sys.sol
```

`MySolutionWrapper` wraps an `AbstractTimeseriesSolution` which already implements the interface.
Since `symbolic_container` will return the wrapped solution, all method calls such as
`is_parameter(sys::MySolutionWrapper, sym)` will be forwarded to `is_parameter(sys.sol, sym)`.

In cases where some methods need to function differently than those of the wrapped type, they can be selectively
defined. For example, suppose `MySolutionWrapper` does not support observed quantities. The following
method can be defined (in addition to the one above):

```julia
is_observed(sys::MySolutionWrapper, sym) = false
```
