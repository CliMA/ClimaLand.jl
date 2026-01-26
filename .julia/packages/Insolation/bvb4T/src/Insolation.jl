module Insolation

using Dates, DelimitedFiles, Interpolations
using Artifacts

include("Parameters.jl")
import .Parameters as IP
const AIP = IP.AbstractInsolationParams

export orbital_params

"""
    OrbitalData

The parameters vary due to Milankovitch cycles.

Orbital parameters from the Laskar 2004 paper are
lazily downloaded from Caltech Box to the
`orbital_parameters_dataset_path(artifact_dir)` path
where `artifact_dir` is the path and filename to save
the artifacts toml file.
"""
struct OrbitalData{E, G, O}
    e_spline_etp::E
    γ_spline_etp::G
    ϖ_spline_etp::O
end

function OrbitalData()
    datapath = joinpath(artifact"laskar2004", "INSOL.LA2004.BTL.csv")
    # We add type annotations here to fix some inference failures.
    Tx = Tuple{Matrix{Float64}, Matrix{AbstractString}}
    x, _ = readdlm(datapath, ',', Float64, header = true)::Tx
    t_range = ((x[1, 1] * 1000):1000:(x[end, 1] * 1000)) # array of every 1 kyr to range of years
    e_spline_etp =
        cubic_spline_interpolation(t_range, x[:, 2]; extrapolation_bc = NaN)
    γ_spline_etp =
        cubic_spline_interpolation(t_range, x[:, 3]; extrapolation_bc = NaN)
    ϖ_spline_etp =
        cubic_spline_interpolation(t_range, x[:, 4]; extrapolation_bc = NaN)

    E = typeof(e_spline_etp)
    G = typeof(γ_spline_etp)
    O = typeof(ϖ_spline_etp)
    return OrbitalData{E, G, O}(e_spline_etp, γ_spline_etp, ϖ_spline_etp)
end

Base.broadcastable(x::OrbitalData) = tuple(x)

e_spline(od, args...) = od.e_spline_etp(args...)
γ_spline(od, args...) = od.γ_spline_etp(args...)
ϖ_spline(od, args...) = od.ϖ_spline_etp(args...)

"""
    orbital_params(od::OrbitalData, dt::FT) where {FT <: Real}

Parameters are interpolated from the values given in the
Laskar 2004 dataset using a cubic spline interpolation.

See [`OrbitalData`](@ref).
"""
function orbital_params(od::OrbitalData, dt::FT) where {FT <: Real}
    return ϖ_spline(od, dt), γ_spline(od, dt), e_spline(od, dt)
end

include("ZenithAngleCalc.jl")
include("InsolationCalc.jl")

# For backwards compatibility with package extensions
if !isdefined(Base, :get_extension)
    include(joinpath("..", "ext", "CreateParametersExt.jl"))
end

end # module
