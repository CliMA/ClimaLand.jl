using Test

import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
using ClimaUtilities.ClimaArtifacts
import Interpolations
import ClimaUtilities
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP

using ClimaLand
import ClimaLand
import ClimaLand.Parameters as LP
const FT = Float64;
regridder_type = :InterpolationsRegridder
extrapolation_bc =
    (Interpolations.Periodic(), Interpolations.Flat(), Interpolations.Flat())
context = ClimaComms.context()
ClimaComms.init(context)

earth_param_set = LP.LandParameters(FT)
radius = FT(6378.1e3)
depth = FT(50)
domain = ClimaLand.Domains.SphericalShell(;
    radius = radius,
    depth = depth,
    nelements = (101, 15),
    npolynomial = 1,
    dz_tuple = FT.((10.0, 0.05)),
)
surface_space = domain.space.surface
subsurface_space = domain.space.subsurface
spatially_varying_soil_params =
    ClimaLand.default_spatially_varying_soil_parameters(
        subsurface_space,
        surface_space,
        FT;
        regridder_type = regridder_type,
        extrapolation_bc = extrapolation_bc,
    )
param_names3d = (:ν, :ν_ss_om, :ν_ss_quartz, :ν_ss_gravel, :K_sat, :S_s, :θ_r)
param_names2d =
    (:PAR_albedo_dry, :NIR_albedo_dry, :PAR_albedo_wet, :NIR_albedo_wet, :f_max)
for p in param_names3d
    @test p ∈ propertynames(spatially_varying_soil_params)
    @test axes(getproperty(spatially_varying_soil_params, p)) ==
          subsurface_space
end

for p in param_names2d
    @test p ∈ propertynames(spatially_varying_soil_params)
    @test axes(getproperty(spatially_varying_soil_params, p)) == surface_space
end
@test :hydrology_cm ∈ propertynames(spatially_varying_soil_params)
hcm = spatially_varying_soil_params.hydrology_cm
@test axes(hcm) == subsurface_space
@test eltype(hcm) == ClimaLand.Soil.vanGenuchten{FT}

clm_parameters = ClimaLand.clm_canopy_parameters(
    surface_space;
    regridder_type = regridder_type,
    extrapolation_bc = extrapolation_bc,
)
param_names = (
    :Ω,
    :rooting_depth,
    :is_c3,
    :Vcmax25,
    :g1,
    :α_PAR_leaf,
    :τ_PAR_leaf,
    :α_NIR_leaf,
    :τ_NIR_leaf,
)
for p in param_names
    @test p ∈ propertynames(clm_parameters)
    @test axes(getproperty(clm_parameters, p)) == surface_space
end

@test :G_Function ∈ propertynames(clm_parameters)
@test axes(getproperty(clm_parameters, :G_Function).χl) == surface_space
# test `use_lowres_clm` on the sphere domain, and then again with a sphere domain with 2x
# the horizontal resolution. Then repeat with plane domains.
@test ClimaLand.use_lowres_clm(surface_space)
domain = ClimaLand.Domains.SphericalShell(;
    radius = radius,
    depth = depth,
    nelements = (202, 2),
    npolynomial = 1,
)
surface_space = domain.space.surface
@test !ClimaLand.use_lowres_clm(surface_space)
domain = ClimaLand.Domains.Plane(;
    longlat = (-117.0, 34.0),
    xlim = (0.0, FT(2e6)),
    ylim = (0.0, 2FT(2e6)),
    nelements = (10, 10),
    npolynomial = 1,
)
surface_space = domain.space.surface
@test ClimaLand.use_lowres_clm(surface_space)
domain = ClimaLand.Domains.Plane(;
    longlat = (-117.0, 34.0),
    xlim = (0.0, FT(2e5)),
    ylim = (0.0, 2FT(2e5)),
    nelements = (10, 10),
    npolynomial = 1,
)
surface_space = domain.space.surface
@test !ClimaLand.use_lowres_clm(surface_space)
