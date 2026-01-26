# Terminology

SymbolicIndexingInterface.jl uses various library-specific terminology throughout its
documentation. This page attempts to provide comprehensive explanations of the terms
used.

## Indexes

An index is an object that defines how to extract specific data from a data structure.
Indexes may be anything from integer indexes into arrays, to custom types that contain
information specific to a particular data structure.

In code samples, an index is typically denoted with the name `idx` or `i`.

## Symbolic variables

Symbolic variables are objects that represent quantities (states, parameters, time, etc.)
or collections of quantities used in a numerical simulation in a more human-accessible
manner. Typically the values of these quantities are stored in custom data structures
and referred to using indexes that do not convey any semantic meaning to users. Symbolic
variables cannot directly be used to access values from these data structures and need
to be translated into indexes.

In code samples, a symbolic variable is typically denoted with the name `sym`.

Symbolic variables are also sometimes referred to as "symbolic indices".

## Index providers

Index providers translate symbolic variables into indexes. In general, an "index" can
be anything from integer indexes into an array, or custom types that define how to
index into complicated data structures. `nothing` is reserved to denote the absence of
an index, in cases where the index provider is unaware of a particular symbolic variable.
It cannot be used as an index. ModelingToolkit.jl systems are examples of index providers.

In code samples, an index provider is typically denoted with the name `indp`.

## Value providers

Value providers store values of symbolic variables. Given an appropriate index from an
index provider, value providers return the value (or values) stored at that index. The
problem, integrator and solution types in SciML are all examples of value providers.
Each value provider is (directly or indirectly) associated with an index provider that
defines the set of valid symbolic variables for that value provider, and the corresponding
indexes.

A value provider may not store values for all symbolic variables in the corresponding index
provider. For example, a parameter object (even a plain `Array` storing parameter values)
is a value provider specifically for the symbolic variables referring to parameters.

In code samples, a value provider is typically denoted with the name `valp`.

!!! note
    It is important to note that an object may be both a value- and index- provider. For
    example, SciML's problem, integrator and solution types are both value- and index-
    providers. This allows for several syntactic improvements. The [`symbolic_container`](@ref)
    function is useful in defining such objects.

### Timeseries objects

Timeseries objects are value providers which implement the [`Timeseries`](@ref) variant of
the [`is_timeseries`](@ref) trait.

### Parameter timeseries objects

Parameter timeseries objects are timeseries objects which implement the
[`Timeseries`](@ref) variant of the [`is_parameter_timeseries`](@ref) trait.
