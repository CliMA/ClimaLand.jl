# ConstructionBase.jl

[`ConstructionBase`](@ref) allows flexible construction and destructuring of objects.
There are two levels of under which this can be done:
### [The raw level](@id the-raw-level)
This is where `Base.fieldnames`, `Base.getfield`, `Base.setfield!` live.
This level is what an object is ultimately composed of including all private details.
At the raw level [`ConstructionBase`](@ref) adds [`constructorof`](@ref) and [`getfields`](@ref).
### [The semantic level](@id the-semantic-level)
This is where `Base.propertynames`, `Base.getproperty` and `Base.setproperty!` live. This level is typically the public interface of a type, it may hide private details and do magic tricks.
At the semantic level [`ConstructionBase`](@ref) adds [`setproperties`](@ref) and [`getproperties`](@ref).


## Interface

```@index
```

```@docs
ConstructionBase
ConstructionBase.constructorof
ConstructionBase.getfields
ConstructionBase.getproperties
ConstructionBase.setproperties
```

## [Tips for designing types](@id type-tips)

When designing types from scratch, it is often possible to structure the types
in such a way that overloading `constructorof` or `setproperties` is unnecessary
in the first place.  It let types in your package work nicely with the ecosystem
built on top of `ConstructionBase` even without explicitly depending on it.
For simple `struct`s whose type parameters can be determined from field values,
`ConstructionBase` works without any customization, provided that the "type-less"
constructor exists.  However, it is often useful or required to have type
parameters that cannot be determined from field values. One way to solve this
problem is to define singleton types that would determine the type parameters:

```jldoctest tips
abstract type OutputBy end
struct Mutating <: OutputBy end
struct Returning <: OutputBy end

struct Add{O <: OutputBy, T}
    outputby::O
    value::T
end

(f::Add{Mutating})(y, x) = y .= x .+ f.value
(f::Add{Returning})(x) = x .+ f.value

add1! = Add(Mutating(), 1)

using ConstructionBase
add2 = constructorof(typeof(add1!))(Returning(), 2)
add2(1)

# output

3
```

`setproperties` works as well:

```jldoctest tips
add3 = setproperties(add2; value=3)
add3(1)

# output

4
```

Note that no overloading of `ConstructionBase` functions was required.
Importantly, this also provides an interface to change type parameters
out-of-the-box:

```jldoctest tips
add3! = setproperties(add3; outputby=Mutating())
add3!([0], 1)

# output

1-element Vector{Int64}:
 4
```

Furthermore, it would work with packages depending on `ConstructionBase` such
as [Setfield.jl](https://github.com/jw3126/Setfield.jl).

```jldoctest tips
using Setfield: @set
add3′ = @set add3!.outputby = Returning()
add3′ === add3

# output

true
```

!!! note

    If it is desirable to keep fields as an implementation detail, combining
    trait functions and
    [`Setfield.FunctionLens`](https://jw3126.github.io/Setfield.jl/latest/#Setfield.FunctionLens)
    may be useful:

    ```jldoctest tips
    OutputBy(x) = typeof(x)
    OutputBy(::Type{<:Add{O}}) where O = O()

    using Setfield: Setfield, @lens
    Setfield.set(add::Add, ::typeof(@lens OutputBy(_)), o::OutputBy) =
        @set add.outputby = o

    obj = (add=add3!,)
    obj′ = @set OutputBy(obj.add) = Returning()
    obj′ === (add=add3,)

    # output

    true
    ```

    ```jldoctest tips
    Setfield.set(::Type{Add{O0, T}}, ::typeof(@lens OutputBy(_)), ::O1) where {O0, T, O1 <: OutputBy} =
        Add{O1, T}

    T1 = typeof(add3!)
    T2 = @set OutputBy(T1) = Returning()
    T2 <: Add{Returning}

    # output

    true
    ```
