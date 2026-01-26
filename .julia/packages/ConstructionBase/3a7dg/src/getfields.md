    getfields(obj) -> NamedTuple
    getfields(obj::Tuple) -> Tuple

Return a `NamedTuple` containing the fields of `obj`. On `Tuples` `getfields` is
the identity function instead, since `Tuple` fields have no symbolic names.

# Examples
```jldoctest
julia> using ConstructionBase

julia> struct S{A,B}
           a::A
           b::B
       end

julia> getfields(S(1,2))
(a = 1, b = 2)

julia> getfields((a=10,b=20))
(a = 10, b = 20)

julia> getfields((4,5,6))
(4, 5, 6)
```

# Specification

`getfields` belongs to the [the raw level](@ref the-raw-level).
Semantically `getfields` boils down to `getfield` and `fieldnames`:
```julia
function getfields(obj::T) where {T}
    fnames = fieldnames(T)
    NamedTuple{fnames}(getfield.(Ref(obj), fnames))
end
```
However the actual implementation can be more optimized. For builtin types, there can also be deviations from this semantics:
* `getfields(::Tuple)::Tuple` since `Tuples` don't have symbolic fieldnames
* There are some types in `Base` that have `undef` fields. Since accessing these results in an error, `getfields` instead just omits these.

# Implementation

The semantics of `getfields` should not be changed for user defined types. It should 
return the raw fields as a `NamedTuple` in the struct order. In other words it should be
equivalent to
```julia
function getfields(obj::T) where {T}
    fnames = fieldnames(T)
    NamedTuple{fnames}(getfield.(Ref(obj), fnames))
end
```
even if that includes private fields of `obj`.
If a change of semantics is desired, consider overloading [`getproperties`](@ref) instead.
See also [`getproperties`](@ref), [`constructorof`](@ref)
