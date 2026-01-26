# IfElse.jl

[![Build Status](https://github.com/SciML/IfElse.jl/workflows/CI/badge.svg)](https://github.com/SciML/IfElse.jl/actions?query=workflow%3ACI)

Sometimes, it's good to have a function form of a conditional. Julia's Base
defines `ifelse` for this, but... psyche, it's not defined in Base but in Core!
While this rarely matters, if you're trying to define a new dispatch for `Core.ifelse`
you will find an interesting error message:

```julia
julia> Base.ifelse(x::Array,y,z) = 2x
ERROR: cannot add methods to a builtin function
Stacktrace:
 [1] top-level scope at none:1
```

Thus no new methods can be added to this function in its current form, blocking
symbolic libraries like [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl)
or vectorization libraries like [LoopVectorization.jl](https://github.com/chriselrod/LoopVectorization.jl)
from extending the operations of this function onto their own types.

To get around this issue, we need a separate version of ifelse that is extendable.
This is difficult to do in the standard case without any inference issues, though
it is being worked on [in this PR](https://github.com/JuliaLang/julia/pull/37343).
In the meantime, to allow for development of `ifelse` dispatches from a common
platform, IfElse.jl has been created to hold this common definition.

Thus IfElse.jl does one thing: it defines `IfElse.ifelse` to call `Core.ifelse`
on its fallbacks, but be free to extend. This way, all of the libraries that
want to extend an `ifelse` function can do so on a common form. If and when
\#37343 is fixed, this function will be deprecated for an extendable `Base.ifelse`,
but for now it will allow for the thriving package ecosystem to continue its
development in this area.
