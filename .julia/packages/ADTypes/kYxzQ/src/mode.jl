"""
    AbstractMode

Abstract supertype for the traits identifying differentiation modes.

# Subtypes

  - [`ForwardMode`](@ref)
  - [`ReverseMode`](@ref)
  - [`ForwardOrReverseMode`](@ref)
  - [`SymbolicMode`](@ref)
"""
abstract type AbstractMode end

"""
    mode(ad::AbstractADType)

Return the differentiation mode of `ad`, as a subtype of [`AbstractMode`](@ref).
"""
function mode end

"""
    ForwardMode

Trait for AD choices that rely on [forward mode](https://en.wikipedia.org/wiki/Automatic_differentiation#Forward_accumulation) algorithmic differentiation or [finite differences](https://en.wikipedia.org/wiki/Finite_difference).

These two paradigms are classified together because they can both efficiently compute Jacobian-vector products.
"""
struct ForwardMode <: AbstractMode end

"""
    ReverseMode

Trait for AD choices that rely on [reverse mode](https://en.wikipedia.org/wiki/Automatic_differentiation#Reverse_accumulation) algorithmic differentiation.
"""
struct ReverseMode <: AbstractMode end

"""
    ForwardOrReverseMode

Trait for AD choices that can work either in [`ForwardMode`](@ref) or [`ReverseMode`](@ref), depending on their configuration.

!!! warning

    This trait should rarely be used, because more precise dispatches to [`ForwardMode`](@ref) or [`ReverseMode`](@ref) should be defined.
"""
struct ForwardOrReverseMode <: AbstractMode end

"""
    SymbolicMode

Trait for AD choices that rely on [symbolic differentiation](https://en.wikipedia.org/wiki/Computer_algebra).
"""
struct SymbolicMode <: AbstractMode end
