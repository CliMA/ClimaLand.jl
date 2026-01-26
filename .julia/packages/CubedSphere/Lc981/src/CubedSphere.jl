module CubedSphere

export sn, cn, conformal_cubed_sphere_mapping, conformal_cubed_sphere_inverse_mapping
export conformal_cubed_sphere_coordinates, optimized_non_uniform_conformal_cubed_sphere_coordinates
export UniformSpacing, GeometricSpacing, ExponentialSpacing

using Printf
using TaylorSeries

struct UniformSpacing end
struct GeometricSpacing end
struct ExponentialSpacing end

include("spherical_geometry.jl")
include("rancic_taylor_coefficients.jl")
include("conformal_cubed_sphere.jl")
include("non_uniform_conformal_mapping_coordinates.jl")

end # module
