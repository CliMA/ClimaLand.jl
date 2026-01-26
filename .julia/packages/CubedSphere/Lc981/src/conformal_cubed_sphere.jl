W_Rancic(Z) = sum(A_Rancic[k] * Z^(k-1) for k in length(A_Rancic):-1:1)
Z_Rancic(W) = sum(B_Rancic[k] * W^(k-1) for k in length(B_Rancic):-1:1)

"""
    conformal_cubed_sphere_mapping(x, y; W_map=W_Rancic)

Conformal mapping from a face of a cube onto the equivalent sector of a sphere with unit radius.

Map the north-pole face of a cube with coordinates ``(x, y)`` onto the equivalent sector of the
sphere with coordinates ``(X, Y, Z)``.

The cube's face oriented normal to ``z``-axis and its coordinates must lie within the
range ``-1 ≤ x ≤ 1``, ``-1 ≤ y ≤ 1`` with its center at ``(x, y) = (0, 0)``. The coordinates
``X, Y`` increase in the same direction as ``x, y``.

The numerical conformal mapping used here is described by [Rancic-etal-1996](@citet).

This is a Julia translation of [MATLAB code from MITgcm](http://wwwcvs.mitgcm.org/viewvc/MITgcm/MITgcm_contrib/high_res_cube/matlab-grid-generator/map_xy2xyz.m?view=markup) that is based on
Fortran 77 code from Jim Purser & Misha Rančić.

Examples
========

The center of the cube's face ``(x, y) = (0, 0)`` is mapped onto ``(X, Y, Z) = (0, 0, 1)``

```jldoctest
julia> using CubedSphere

julia> conformal_cubed_sphere_mapping(0, 0)
(0, 0, 1.0)
```

and the edge of the cube's face at ``(x, y) = (1, 1)`` is mapped onto ``(X, Y, Z) = (√3/3, √3/3, √3/3)``

```jldoctest
julia> using CubedSphere

julia> conformal_cubed_sphere_mapping(1, 1)
(0.5773502691896256, 0.5773502691896256, 0.5773502691896257)
```

# References

* [Rancic-etal-1996](@cite) Rančić et al., *Q. J. R. Meteorol.*, (1996).
"""
function conformal_cubed_sphere_mapping(x, y; W_map=W_Rancic)
    (abs(x) > 1 || abs(y) > 1) && throw(ArgumentError("(x, y) must lie within [-1, 1] x [-1, 1]"))

    X = xᶜ = abs(x)
    Y = yᶜ = abs(y)

    kxy = yᶜ > xᶜ

    xᶜ = 1 - xᶜ
    yᶜ = 1 - yᶜ

    kxy && (xᶜ = 1 - Y)
    kxy && (yᶜ = 1 - X)

    Z = ((xᶜ + im * yᶜ) / 2)^4
    W = W_map(Z)

    im³ = im^(1/3)
    ra = √3 - 1
    cb = -1 + im
    cc = ra * cb / 2

    W = im³ * (W * im)^(1/3)
    W = (W - ra) / (cb + cc * W)
    X, Y = reim(W)

    H = 2 / (1 + X^2 + Y^2)
    X = X * H
    Y = Y * H
    Z = H - 1

    if kxy
        X, Y = Y, X
    end

    y < 0 && (Y = -Y)
    x < 0 && (X = -X)

    # Fix truncation for x = 0 or y = 0.
    x == 0 && (X = 0)
    y == 0 && (Y = 0)

    return X, Y, Z
end

"""
    conformal_cubed_sphere_inverse_mapping(X, Y, Z; Z_map=Z_Rancic)

Inverse mapping for conformal cube sphere for quadrant of North-pole face in which `X` and `Y` are both
positive. All other mappings to other cube face coordinates can be recovered from rotations of this map.
There a 3 other quadrants for the north-pole face and five other faces for a total of twenty-four quadrants.
Because of symmetry only the reverse for a single quadrant is needed. Because of branch cuts and the complex
transform the inverse mappings are multi-valued in general, using a single quadrant case allows a simple
set of rules to be applied.

The mapping is valid for the cube face quadrant defined by ``0 < x < 1`` and ``0 < y < 1``, where a full cube
face has extent ``-1 < x < 1`` and ``-1 < y < 1``. The quadrant for the mapping is from a cube face that has
"north-pole" at its center ``(x=0, y=0)``. i.e., has `X, Y, Z = (0, 0, 1)` at its center. The valid ranges of `X`
and `Y` for this mapping and convention are a quadrant defined be geodesics that connect the points A, B, C and D,
on the shell of a sphere of radius ``R`` with `X`, `Y` coordinates as follows

```
A = (0, 0)
B = (√2, 0)
C = (√3/3, √3/3)
D = (0, √2)
```
"""
function conformal_cubed_sphere_inverse_mapping(X, Y, Z; Z_map=Z_Rancic)
    H  = Z + 1
    Xˢ = X / H
    Yˢ = Y / H
    ω  = Xˢ + im * Yˢ

    ra = √3 - 1
    cb = -1 + im
    cc = ra * cb / 2
    ω⁰ = (ω * cb + ra) / (1 - ω * cc)
    W⁰ = im * ω⁰^3 * im
    Z  = Z_map(W⁰)
    z  = 2 * Z^(1/4)
    x, y = reim(z)

    kxy = abs(y) > abs(x)
    xx = x
    yy = y
    !kxy && ( x = 1 - abs(yy) )
    !kxy && ( y = 1 - abs(xx) )

    xf = x
    yf = y

    ( X < Y ) && ( xf = y )
    ( X < Y ) && ( yf = x )

    x = xf
    y = yf

    return x, y
end

"""
    cube_to_sphere(x, y)

Maps the coordinates ``(x, y) ∈ [-1, 1] × [-1, 1]`` of the face of a cube to the
3D coordinates ``X, Y, Z`` on the sphere via [`conformal_cubed_sphere_mapping`](@ref),
that is:

    X[i, j], Y[i, j], Z[i, j] = conformal_cubed_sphere_mapping(x[i], y[j])
"""
function cube_to_sphere(x, y)
    X = zeros(length(x), length(y))
    Y = zeros(length(x), length(y))
    Z = zeros(length(x), length(y))

    for (j, y′) in enumerate(y), (i, x′) in enumerate(x)
        X[i, j], Y[i, j], Z[i, j] = conformal_cubed_sphere_mapping(x′, y′)
    end

    return X, Y, Z
end
