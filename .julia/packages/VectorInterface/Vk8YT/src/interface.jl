# Vector Interface
#------------------
# each of the following methods need to be implemented by custom types

"""
    scalartype(x)

Returns the type of scalar over which the vector-like object `x` behaves as a vector, e.g.
the type of scalars with which `x` could be scaled in-place. This function should also work
in the type domain, i.e. if `x` is a vector like object, then also `scalartype(typeof(x))`
should work.

!!! note
    New types should only implement the method with the argument in the type domain.
"""
function scalartype end
scalartype(x) = scalartype(typeof(x))
# implementation should be in type domain

## zerovector
"""
    zerovector(x, [S::Type{<:Number} = scalartype(x)])

Returns a zero vector in the vector space of `x`. Optionally, a modified scalar type `S` for
the resulting zero vector can be specified.

!!! note
    New types should only implement the two-argument version, if applicable.

Also see: [`zerovector!`](@ref), [`zerovector!!`](@ref)
"""
function zerovector end
zerovector(x) = zerovector(x, scalartype(x))

"""
    zerovector!(x)

Modifies `x` in-place to become the zero vector.

Also see: [`zerovector`](@ref), [`zerovector!!`](@ref)
"""
function zerovector! end

"""
    zerovector!!(x, [S::Type{<:Number} = scalartype(x)])

Construct a zero vector in the vector space of `x`, thereby trying to overwrite and thus
recycle `x` when possible. Optionally, a modified scalar type `S` for the resulting zero
vector can be specified.

!!! note
    New types should only implement the one-argument version. The two-argument version
    amounts to `S == scalartype(x) ? zerovector!!(x) : zerovector(x, S)`

Also see: [`zerovector`](@ref), [`zerovector!`](@ref)
"""
function zerovector!! end
zerovector!!(x, S::Type{<:Number}) = S == scalartype(x) ? zerovector!!(x) : zerovector(x, S)

"""
    scale(x, α::Number)

Computes the new vector-like object obtained from scaling `x` with the scalar `α`.

Also see: [`scale!`](@ref) and [`scale!!`](@ref)
"""
function scale end

"""
    scale!(x, α::Number) -> x
    scale!(y, x, α::Number) -> y

Rescale `x` with the scalar coefficient `α`, thereby overwrite and thus recylcing the
contents of `x` (in the first form) or `y` (in the second form). This is only possible if
`x`, respectively `y` is mutable, and if the scalar types involed are compatible and can
be promoted and converted.

Also see: [`scale`](@ref) and [`scale!!`](@ref)
"""
function scale! end

"""
    scale!!(x, α::Number)
    scale!!(y, x, α::Number)

Rescale `x` with the scalar coefficient `α`, thereby trying to overwrite and thus recylce
the contents of `x` (in the first form) or `y` (in the second form).  A new object will be
created when this fails due to immutability, incompatible scalar types or
incommensurate sizes.

Also see: [`scale`](@ref) and [`scale!`](@ref)
"""
function scale!! end

"""
    add(y, x, [α::Number = 1, β::Number = 1])

Add `y` and `x`, or more generally construct the linear combination `y * β + x * α`.

!!! note
    New types should only implement the four-argument version. When desired, the two- and
    three-argument version can be distinguished by specializing on `α` and `β` being of
    type `One`.

See also: [`add!`](@ref) and [`add!!`](@ref)
"""
add(y, x) = add(y, x, One(), One())
add(y, x, α::Number) = add(y, x, α, One())

"""
    add!(y, x, [α::Number = One(), β::Number = One()])

Add `y` and `x`, or more generally construct the linear combination `y * β + x * α`, storing
the result in `y`. This will error in case of incompatible scalar types or incommensurate sizes.

!!! note
    New types should only implement the four-argument version. When desired, the two- and
    three-argument version can be distinguished by specializing on `α` and `β` being of
    type `One`.

See also: [`add`](@ref) and [`add!!`](@ref)
"""
add!(y, x) = add!(y, x, One(), One())
add!(y, x, α::Number) = add!(y, x, α, One())

"""
    add!!(y, x, [α::Number = 1, β::Number = 1])

Add `y` and `x`, or more generally construct the linear combination `y * β + x * α`, thereby
trying to store the result in `y`. A new object will be created when this fails due to
immutability, incompatible scalar types or incommensurate sizes.

!!! note
    New types should only implement the four-argument version. When desired, the two- and
    three-argument version can be distinguished by specializing on `α` and `β` being of
    type `One`.

See also: [`add`](@ref) and [`add!`](@ref)
"""
add!!(y, x) = add!!(y, x, One(), One())
add!!(y, x, α::Number) = add!!(y, x, α, One())

# Inner product
"""
    inner(x, y)

Compute the inner product between `x` and `y`.
"""
function inner end
