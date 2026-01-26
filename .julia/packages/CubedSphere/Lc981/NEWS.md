CubedSphere.jl Release Notes
===============================

v0.3.3
------

Introduces conformal cubed-sphere coordinates which originate from a non-uniformly
spaced horizontal grid (ξ, η) ∈ [-1, 1] x [-1, 1].


v0.3.0
------

We moved `conformal_map_tylor_coefficients.jl` and `complex_jacobi_ellptic.jl`
from main package to `sandbox`. This allowed us to remove most dependencies of
`CubedSphere.jl`. `sandbox` now contains a Julia environment with the packages
needed to execute the scripts in that folder.

