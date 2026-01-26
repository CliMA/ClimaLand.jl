# The SciMLStructure Interface

## Core Interface Definitions

### `isscimlstructure` Definition

```julia
isscimlstructure(p)::Bool
ismutablescimlstructure(p)::Bool
```

Returns whether the object satisfies the SciMLStructure interface. Defaults to `false` and types
are required to opt-into the interface.

### `canonicalize` Definition

```julia
canonicalize(::AbstractPortion, p::T1) -> values::T2, repack, aliases::Bool
repack(new_values::T2) -> p::T1 # with values replaced with new_values
replace(::AbstractPortion, p::T1, new_values) -> p::T1
replace!(::AbstractPortion, p::T1, new_values)::Nothing # Requires mutable
```

### Portion Definitions

The core function of the interface is the `canonicalize` function. `canonicalize` allows the user to define
to the solver how to represent the given "portion" in a standard `AbstractVector` type which allows for
interfacing with standard tools like linear algebra in an efficient manner. The type of portions which
are defined are:

  - Tunable: the tunable values/parameters, i.e. the values of the structure which are supposed to be considered
    non-constant when used in the context of an inverse problem solve. For example, this is the set of
    parameters to be optimized during a parameter estimation of an ODE.
    
      + Tunable parameters are expected to return an `AbstractVector` of unitless values.
      + Tunable parameters are expected to be constant during the solution of the ODE.

  - Constants: the values which are to be considered constant by the solver, i.e. values which are not estimated
    in an inverse problem and which are unchanged in any operation by the user as part of the solver's usage.
  - Caches: the stored cache values of the struct, i.e. the values of the structure which are used as intermediates
    within other computations in order to avoid extra allocations.
  - Discrete: the discrete portions of the state captured inside of the structure. For example, discrete values
    stored outside of the `u` in the parameters to be modified in the callbacks of an ODE.
    
      + Any parameter that is modified inside of callbacks should be considered Discrete.

## Definitions for Base Objects

  - `Vector`: returns an aliased version of itself as `Tunable`, and an empty vector matching type for `Constants`,
    `Caches`, and `Discrete`.
  - `Array`: returns the `vec(p)` aliased version of itself as `Tunable`, and an empty vector matching type for `Constants`,
    `Caches`, and `Discrete`.
