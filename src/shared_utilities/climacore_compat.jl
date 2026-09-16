# ClimaCore 1.0.0 dropped the `Spaces.issubspace` methods that let a field on a
# level of an extruded space broadcast against a field on the horizontal space
# they share, which ClimaLand relies on for surface boundary conditions. The
# `Spaces.horizontal_grid`/`Spaces.vertical_grid` docstrings still describe them,
# so restore them here for any ClimaCore version where they are absent.

_issubspace_is_generic(::Type{T}) where {T} =
    which(Spaces.issubspace, Tuple{T, T}) ===
    which(Spaces.issubspace, Tuple{Spaces.AbstractSpace, Spaces.AbstractSpace})

if _issubspace_is_generic(Spaces.AbstractSpectralElementSpace)
    # Slab spaces carry their local geometry directly and have no grid.
    Spaces.issubspace(
        space1::Spaces.SpectralElementSpaceSlab,
        space2::Spaces.SpectralElementSpaceSlab,
    ) = space1 == space2
    Spaces.issubspace(
        space1::Spaces.AbstractSpectralElementSpace,
        space2::Spaces.AbstractSpectralElementSpace,
    ) =
        Spaces.horizontal_grid(Spaces.grid(space1)) ===
        Spaces.horizontal_grid(Spaces.grid(space2))
end

if _issubspace_is_generic(Spaces.FiniteDifferenceSpace)
    Spaces.issubspace(
        space1::Spaces.FiniteDifferenceSpace,
        space2::Spaces.FiniteDifferenceSpace,
    ) =
        Spaces.vertical_grid(Spaces.grid(space1)) ===
        Spaces.vertical_grid(Spaces.grid(space2))
end

if _issubspace_is_generic(Spaces.MultiPointSpace)
    Spaces.issubspace(
        space1::Spaces.MultiPointSpace,
        space2::Spaces.MultiPointSpace,
    ) =
        Spaces.horizontal_grid(Spaces.grid(space1)) ===
        Spaces.horizontal_grid(Spaces.grid(space2))
end
