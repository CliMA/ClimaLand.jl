module CommonSolve

"""
```julia
CommonSolve.solve(args...; kwargs...)
```

Solves an equation or other mathematical problem using the algorithm
specified in the arguments. Generally, the interface is:

```julia
CommonSolve.solve(prob::ProblemType,alg::SolverType; kwargs...)::SolutionType
```

where the keyword arguments are uniform across all choices of algorithms.

By default, `solve` defaults to using `solve!` on the iterator form, i.e.:

```julia
solve(args...; kwargs...) = solve!(init(args...; kwargs...))
```
"""
solve(args...; kwargs...) = solve!(init(args...; kwargs...))

"""
```julia
CommonSolve.solve!(iter)
```

Solves an equation or other mathematical problem using the algorithm
specified in the arguments. Generally, the interface is:

```julia
iter = CommonSolve.init(prob::ProblemType,alg::SolverType; kwargs...)::IterType
CommonSolve.solve!(iter)::SolutionType
```

where the keyword arguments are uniform across all choices of algorithms.
The `iter` type will be different for the different problem types.
"""
function solve! end

"""
```julia
iter = CommonSolve.init(args...; kwargs...)
```

Solves an equation or other mathematical problem using the algorithm
specified in the arguments. Generally, the interface is:

```julia
iter = CommonSolve.init(prob::ProblemType,alg::SolverType; kwargs...)::IterType
CommonSolve.solve!(iter)::SolutionType
```

where the keyword arguments are uniform across all choices of algorithms.
The `iter` type will be different for the different problem types.
"""
function init end

"""
```julia
CommonSolve.step!(iter, args...; kwargs...)
```

Progress the iterator object (the one returned by `CommonSolve.init`).
The additional arguments typically describe how much to progress
the iterator for, and are implementation-specific.
"""
function step! end

end # module
