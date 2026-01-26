module SphericalGeometry

export lat_lon_to_x, lat_lon_to_y, lat_lon_to_z, lat_lon_to_cartesian, cartesian_to_latitude, cartesian_to_longitude,
    cartesian_to_lat_lon, spherical_distance, spherical_area_triangle, spherical_area_quadrilateral,
    spherical_quadrilateral_vertices, compute_cell_areas

using Distances
using LinearAlgebra

"""
    cartesian_to_lat_lon(x, y, z)
    cartesian_to_lat_lon(X)

Convert 3D Cartesian coordinates `(x, y, z)` or a 3-element Cartesian vector `X = (x, y, z)` on the sphere to
latitude–longitude. Returns a tuple `(latitude, longitude)` in degrees.

- Latitude is the angle measured from the equatorial plane (`z = 0`).
- Longitude is measured anti-clockwise (eastward) from the `x`-axis (`y = 0`) about the `z`-axis.

# Arguments
- `x, y, z`: Cartesian coordinates (numbers), **or**
- `X`: 3-element Cartesian vector.

# Returns
- `(latitude, longitude)`: Latitude and longitude angles in degrees.

# Examples
Find latitude–longitude of the North Pole:

```jldoctest 1
julia> using CubedSphere.SphericalGeometry

julia> x, y, z = (0, 0, 6.4e6);  # Cartesian coordinates of the North Pole [in meters]

julia> cartesian_to_lat_lon(x, y, z)
(90.0, 0.0)
```
Let's confirm that for few points on the unit sphere we get the answers we expect.

```jldoctest 1
julia> cartesian_to_lat_lon(√2/4, -√2/4, √3/2)
(59.99999999999999, -45.0)

julia> cartesian_to_lat_lon(-√6/4, √2/4, -√2/2)
(-45.0, 150.0)
```
"""
cartesian_to_lat_lon(x, y, z) = cartesian_to_latitude(x, y, z), cartesian_to_longitude(x, y, z)

function cartesian_to_lat_lon(X)
    x, y, z = X
    return cartesian_to_lat_lon(x, y, z)
end

"""
    cartesian_to_latitude(x, y, z)

Convert Cartesian coordinates `(x, y, z)` to latitude (in degrees) on the sphere.
"""
cartesian_to_latitude(x, y, z) = atand(z, sqrt(x*x + y*y))

"""
    cartesian_to_longitude(x, y, z)

Convert Cartesian coordinates `(x, y, z)` to longitude (in degrees) on the sphere.
"""
cartesian_to_longitude(x, y, z) = atand(y, x)

"""
lat_lon_to_cartesian(φ, λ; radius = 1, check_latitude_bounds = true)

Convert `(latitude, longitude)` coordinates (in degrees) to Cartesian coordinates `(x, y, z)` on the sphere.

# Arguments
- `φ`: Latitude in degrees.
- `λ`: Longitude in degrees.
- `radius`: Sphere radius (optional). Default is `1`.
- `check_latitude_bounds`: If `true`, raises an error when `|φ| > 90`. Set to `false` to disable this validation.

# Returns
- `(x, y, z)`: Cartesian coordinates on the sphere.

# Examples
Find the Cartesian coordinates of the North Pole on a unit sphere:

```jldoctest 1
julia> using CubedSphere.SphericalGeometry

julia> lat_lon_to_cartesian(90, 0)
(0.0, 0.0, 1.0)
```
Find the Cartesian coordinates of a point on the equator with longitude 90°E:

```jldoctest 1
julia> lat_lon_to_cartesian(0, 90)
(0.0, 1.0, 0.0)
```
"""
function lat_lon_to_cartesian(φ, λ; radius = 1, check_latitude_bounds = true)
    check_latitude_bounds && abs(φ) > 90 && error("Latitude φ must be within -90 ≤ φ ≤ 90 degrees.")
    return (lat_lon_to_x(φ, λ; radius), lat_lon_to_y(φ, λ; radius), lat_lon_to_z(φ; radius))
end

"""
lat_lon_to_x(φ, λ; radius = 1)

Convert (latitude, longitude) coordinates (in degrees) to Cartesian coordinate x on the sphere.
"""
lat_lon_to_x(φ, λ; radius = 1) = radius * cosd(λ) * cosd(φ)

"""
lat_lon_to_y(φ, λ; radius = 1)

Convert (latitude, longitude) coordinates (in degrees) to Cartesian coordinate y on the sphere.
"""
lat_lon_to_y(φ, λ; radius = 1) = radius * sind(λ) * cosd(φ)

"""
lat_lon_to_z(φ; radius = 1)

Convert (latitude, longitude) coordinates (in degrees) to Cartesian coordinate z on the sphere.
"""
lat_lon_to_z(φ; radius = 1) = radius * sind(φ)

"""
    spherical_distance(a₁, a₂; radius = 1)

Compute the great-circle distance between two Cartesian points `a₁` and `a₂` on a sphere of `radius`.

For `radius = 1`, this is equivalent to the central angle (in radians) between the two points.

# Arguments
- `a₁`, `a₂`: 3-element Cartesian vectors on the sphere.
- `radius`: Sphere radius (optional). Default is `1`.

# Returns
- Great-circle distance between `a₁` and `a₂`

# Notes
- Both `a₁` and `a₂` must lie on the surface of the same sphere (i.e., have the same norm).

# Examples
```jldoctest 1
julia> using CubedSphere.SphericalGeometry

julia> a₁ = (1.0, 0.0, 0.0);  # point on unit sphere
       a₂ = (0.0, 1.0, 0.0);  # 90° away

julia> spherical_distance(a₁, a₂)
1.5707963267948968
```
"""
function spherical_distance(a₁, a₂; radius=1)
    (sum(a₁.^2) ≈ sum(a₂.^2)) || error("a₁ and a₂ must have same norm")

    φ₁, λ₁ = cartesian_to_lat_lon(a₁)
    φ₂, λ₂ = cartesian_to_lat_lon(a₂)

    return haversine((λ₁, φ₁), (λ₂, φ₂), radius)
end

"""
    spherical_area_triangle(a::Number, b::Number, c::Number; radius = 1)

Compute the area of a spherical triangle on a sphere of `radius`, given its three side
lengths `a`, `b`, and `c` (in radians).

For a unit sphere (`radius = 1`), the area equals the spherical excess `E = A + B + C − π`,
where `A`, `B`, and `C` are the triangle’s interior angles. For a sphere of radius `R`, the
area is `R² * E`.

# Arguments
- `a`, `b`, `c`: Side lengths of the spherical triangle, measured as central angles (in radians).
- `radius`: Sphere radius (optional). Default is `1`.

# Returns
- The area of the spherical triangle on a sphere of `radius`.

# Notes
- Euler (1778) and Lagrange (1798) showed that the spherical excess `E` on the unit sphere is computed as

  ```math
  \\tan\\frac{E}{2} =
  \\frac{\\sqrt{1 - \\cos^2 a - \\cos^2 b - \\cos^2 c + 2 \\cos a \\cos b \\cos c}}{1 + \\cos a + \\cos b + \\cos c}.
  ```

References
==========

* Euler, L. (1778) De mensura angulorum solidorum, Opera omnia, 26, 204-233 (Orig. in Acta adac. sc. Petrop. 1778)
* Lagrange,  J.-L. (1798) Solutions de quelques problèmes relatifs au triangles sphériques, Oeuvres, 7, 331-359.

# Examples
```jldoctest 1
julia> using CubedSphere.SphericalGeometry

julia> a = b = c = π/2;  # Right spherical triangle with 90° sides on unit sphere

julia> spherical_area_triangle(a, b, c)
1.5707963267948966

julia> spherical_area_triangle(a, b, c; radius = 6371e3)
6.375805898872353e13
```
"""
function spherical_area_triangle(a::Number, b::Number, c::Number; radius=1)
    cosa, cosb, cosc = cos(a), cos(b), cos(c)

    tan½E = sqrt(1 - cosa^2 - cosb^2 - cosc^2 + 2cosa * cosb * cosc) / (1 + cosa + cosb + cosc)

    E_unit = 2atan(tan½E)       # area on unit sphere

    return (radius^2) * E_unit  # physical area
end

"""
    spherical_area_triangle(a₁, a₂, a₃; radius = 1)

Compute the area of a spherical triangle on a sphere of `radius`, given its three vertex position vectors
`a₁`, `a₂`, and `a₃` in 3D Cartesian coordinates. The origin is assumed to be at the center of the sphere.

For a unit sphere (`radius = 1`), the area equals the spherical excess `E`. For a sphere of radius `R`,
the area is `R² * E`.

# Arguments
- `a₁`, `a₂`, `a₃`: 3-element Cartesian vectors representing the vertices of the spherical triangle.
  All three must lie on the same sphere.
- `radius`: Sphere radius (optional). Default is `1`.

# Returns
- The area of the spherical triangle on a sphere of `radius`.

# Notes
- This function generalizes the classical Euler–Lagrange formula for spherical excess by expressing the quantity
  ```math
  P = \\sqrt{1 - \\cos^2 a - \\cos^2 b - \\cos^2 c + 2 \\cos a \\cos b \\cos c}
  ```
  in terms of the scalar triple product
  ```math
  P = |𝐚₁ ⋅ (𝐚₂ × 𝐚₃)|
  ```
  where `a`, `b`, and `c` are the side lengths of the spherical triangle formed by the vertices `a₁`, `a₂`, and `a₃`.
  The above formula was first derived by Eriksson (1990).

- The inputs `a₁`, `a₂`, and `a₃` must have the same norm.

References
==========

* Euler, L. (1778) De mensura angulorum solidorum, Opera omnia, 26, 204-233 (Orig. in Acta adac. sc. Petrop. 1778)
* Lagrange,  J.-L. (1798) Solutions de quelques problèmes relatifs au triangles sphériques, Oeuvres, 7, 331-359.
* Eriksson, F. (1990) On the measure of solid angles, Mathematics Magazine, 63 (3), 184-187,
doi:10.1080/0025570X.1990.11977515

# Examples
```jldoctest 1
julia> using CubedSphere.SphericalGeometry

julia> a₁ = [1.0, 0.0, 0.0];
       a₂ = [0.0, 1.0, 0.0];
       a₃ = [0.0, 0.0, 1.0];

julia> spherical_area_triangle(a₁, a₂, a₃)
1.5707963267948966

julia> spherical_area_triangle(a₁ .* 6.371e6, a₂ .* 6.371e6, a₃ .* 6.371e6; radius = 6.371e6)
6.375805898872353e13
```
"""
function spherical_area_triangle(a₁, a₂, a₃; radius=1)
    a₁, a₂, a₃ = collect(a₁), collect(a₂), collect(a₃)
    r1, r2, r3 = sqrt(sum(a₁.^2)), sqrt(sum(a₂.^2)), sqrt(sum(a₃.^2))
    (r1 ≈ r2 && r2 ≈ r3) || error("a₁, a₂, a₃ must lie on the same sphere")

    # Use only directions to compute unit-sphere area (scale later)
    û₁, û₂, û₃ = a₁/r1, a₂/r2, a₃/r3
    tan½E = abs(dot(û₁, cross(û₂, û₃))) / (1 + dot(û₁, û₂) + dot(û₂, û₃) + dot(û₁, û₃))
    E_unit = 2atan(tan½E)

    return (radius^2) * E_unit
end

"""
    spherical_area_quadrilateral(a₁, a₂, a₃, a₄; radius = 1)

Compute the area of a spherical quadrilateral on a sphere of `radius`, given the position of
its four vertices as vectors `a₁`, `a₂`, `a₃`, and `a₄` in 3D Cartesian coordinates.
The origin is assumed to be at the center of the sphere.

The quadrilateral area is evaluated as half the sum of the areas of all four spherical triangles
formed by the vertices. This approach avoids the need to explicitly choose a diagonal that splits
the quadrilateral into two non-overlapping triangles.

# Arguments
- `a₁`, `a₂`, `a₃`, `a₄`: 3-element Cartesian vectors representing the four vertices of the spherical quadrilateral.
  All four must lie on the same sphere.
- `radius`: Sphere radius (optional). Default is `1`.

# Returns
- The area of the spherical quadrilateral on a sphere of `radius`.

# Note
- This method is numerically robust and works for convex spherical quadrilaterals without requiring explicit diagonal
  selection.

# Examples
```jldoctest 1
julia> using CubedSphere.SphericalGeometry

julia> a₁ = [1.0, 0.0, 0.0];
       a₂ = [0.0, 1.0, 0.0];
       a₃ = [0.0, 0.0, 1.0];
       a₄ = [1.0, 1.0, 0.0] ./ √2;  # mid-edge on unit sphere

julia> spherical_area_quadrilateral(a₁, a₂, a₃, a₄)
1.5707963267948966

julia> R = 6.371e6;

julia> spherical_area_quadrilateral(a₁ .* R, a₂ .* R, a₃ .* R, a₄ .* R; radius = R)
6.375805898872353e13
```
"""
spherical_area_quadrilateral(a₁, a₂, a₃, a₄; radius=1) =
    0.5 * (spherical_area_triangle(a₁, a₂, a₃; radius) +
           spherical_area_triangle(a₁, a₂, a₄; radius) +
           spherical_area_triangle(a₁, a₃, a₄; radius) +
           spherical_area_triangle(a₂, a₃, a₄; radius))

"""
    spherical_quadrilateral_vertices(X, Y, Z, i, j)

Returns the four Cartesian vertex vectors of the spherical grid cell whose corners are indexed by `(i, j)`, `(i+1, j)`,
`(i+1, j+1)`, and `(i, j+1)` in the arrays `X`, `Y`, and `Z`. Each of `X`, `Y`, and `Z` is a 2D array of size `(Nx, Ny)`
holding the Cartesian coordinates of grid vertices on the sphere, such that the point at `(i, j)` is
`(X[i, j], Y[i, j], Z[i, j])`.

# Arguments
- `X`, `Y`, `Z`: `(Nx, Ny)` arrays of Cartesian coordinates on the sphere.
- `i`, `j`: Indices of the lower-left corner of the cell (in array order).

# Returns
- `(a₁, a₂, a₃, a₄)`: The four 3-element Cartesian vertex vectors at `(i, j)`, `(i+1, j)`, `(i+1, j+1)`, and `(i, j+1)`.
"""
function spherical_quadrilateral_vertices(X, Y, Z, i, j)
    x₁ = X[i, j]
    y₁ = Y[i, j]
    z₁ = Z[i, j]
    a₁ = [x₁, y₁, z₁]
    x₂ = X[i+1, j]
    y₂ = Y[i+1, j]
    z₂ = Z[i+1, j]
    a₂ = [x₂, y₂, z₂]
    x₃ = X[i+1, j+1]
    y₃ = Y[i+1, j+1]
    z₃ = Z[i+1, j+1]
    a₃ = [x₃, y₃, z₃]
    x₄ = X[i, j+1]
    y₄ = Y[i, j+1]
    z₄ = Z[i, j+1]
    a₄ = [x₄, y₄, z₄]

    return a₁, a₂, a₃, a₄
end

"""
    compute_cell_areas(X, Y, Z; radius = 1)

Compute the spherical surface areas of all quadrilateral cells in a spherical grid (e.g., a conformal cubed-sphere
panel), defined by the Cartesian coordinate arrays `X`, `Y`, and `Z`. Each of `X`, `Y`, and `Z` is a 2D array of size
`(Nx, Ny)` holding the Cartesian coordinates of the grid vertices on the sphere, such that the point at `(i, j)` is
`(X[i, j], Y[i, j], Z[i, j])`. The grid therefore contains `(Nx−1) × (Ny−1)` spherical quadrilateral cells.

For each quadrilateral cell with vertices `(i, j)`, `(i+1, j)`, `(i+1, j+1)`, and `(i, j+1)`, the function computes the
spherical area using `spherical_area_quadrilateral` and stores the result in a 2D array. If `radius = 1`, the returned
areas correspond to the unit sphere. For `radius ≠ 1`, the physical areas are returned.

# Arguments
- `X`, `Y`, `Z`: `(Nx, Ny)` arrays of Cartesian coordinates of grid vertices on the sphere.
- `radius`: Sphere radius (optional). Default is `1`.

# Returns
- `cell_areas`: An `(Nx−1, Ny−1)` array of quadrilateral cell areas.
  - For `radius = 1`, these are unit-sphere areas.
  - For `radius ≠ 1`, these are physical areas (e.g., in m²).

# Notes
- The function assumes that `X`, `Y`, and `Z` represent vertices lying on the same sphere.

# Examples
```jldoctest 1
julia> using CubedSphere.SphericalGeometry

julia> Nx, Ny = 3, 3;

julia> lons = range(-π/6, π/6, length = Nx);

julia> lats = range(-π/6,  π/6, length = Ny);

julia> X = [cos(φ)*cos(λ) for λ in lons, φ in lats];

julia> Y = [cos(φ)*sin(λ) for λ in lons, φ in lats];

julia> Z = [sin(φ)        for λ in lons, φ in lats];

julia> A = compute_cell_areas(X, Y, Z);

julia> size(A)
(2, 2)

julia> all(isapprox.(A, fill(A[1, 1], 2, 2); rtol=1e-12))
true

julia> isapprox(A[1, 1], 0.26636308214247195; rtol=1e-12)
true

julia> R = 6.371e6;

julia> A_R = compute_cell_areas(X .* R, Y .* R, Z .* R; radius = R);

julia> isapprox.(A_R, A .* R^2; rtol=1e-9) |> all
true
```
"""
function compute_cell_areas(X, Y, Z; radius=1)
    Nx, Ny = size(X)
    cell_areas = zeros(Nx-1, Ny-1)

    for j in 1:Ny-1, i in 1:Nx-1
        a₁, a₂, a₃, a₄ = spherical_quadrilateral_vertices(X, Y, Z, i, j)
        cell_areas[i, j] = spherical_area_quadrilateral(a₁, a₂, a₃, a₄; radius)
    end

    return cell_areas
end

end # module
