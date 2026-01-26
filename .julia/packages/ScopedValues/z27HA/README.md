# ScopedValues.jl

Implement dynamically scoped values for Julia.
This package was primarly inspired by [ContextVariablesX.jl](https://github.com/tkf/ContextVariablesX.jl),
it's corresponding [pull-request to Julia base](https://github.com/JuliaLang/julia/pull/35833), and
Java's [JEPS446](https://openjdk.org/jeps/446).

It has been submitted to Julia as a [JULEP](https://github.com/JuliaLang/julia/pull/50958).

# Usage

```julia
using ScopedValues

const svar = ScopedValue(1)

@show svar[]

# Enter a new dynamic scope and set value
with(svar => 2) do
    @show svar
end

# While a ScopedValue contents are constant
# You can store mutable values.
const svar_dict = ScopedValue(Dict())

# Important we are using `merge` to "unshare" the mutable values
# across the different views of the same scoped value.
with(svar_dict => merge(svar_dict[], Dict(:a => 10))) do
    @show svar_dict[][:a]
end
```

scoped values are inherited across tasks:

```julia
const LEVEL = ScopedValue(:GUEST)

function serve(request, response)
    level = isAdmin(request) ? :ADMIN : :GUEST
    with(LEVEL => level) do
        Threads.@spawn handle(request, respone)
    end
end

function open(connection::Database)
    level = LEVEL[]
    if level !== :ADMIN
        error("Access disallowed")
    end
    # ... open connection
end

function handle(request, response)
    open(Database(#=...=#))
    # ...
end
```

## Usage with logging

Before Julia 1.11, this package implements scoped variables by using a special "payload" logger
(just like [ContextVariablesX.jl](https://github.com/tkf/ContextVariablesX.jl)). Therefore,
changing the logger from within a dynamic scope can clobber the scoped values. To avoid this, use

```julia
ScopedValues.with_logger
```

as a replacement for `Logging.with_logger`. Likewise, if you need to access the logger, use

```julia
ScopedValues.current_logger()
```

which will unwrap the special ScopedValues's "payload" logger to provide you with the actual logging-logger.

Note that there is no issue with setting a logger *outside* of the dynamic scope. For example, a common usage is to set the global logger at the start of the program (and to not change it mid-execution). In this case, there is no issue clash with ScopedValues.jl, and e.g. `Logging.global_logger()` may be used to set the logger. There is only concern when modifying or accessing the logger from within a dynamic scope provided by ScopedValues.jl (i.e. within a `with` block).

There is also no issue with logging to a logger (e.g. `@info`, etc), within a dynamic scope or elsewhere. ScopedValues "payload" logger will pass log messages through to the active logger.
