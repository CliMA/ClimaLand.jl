using Test
import ClimaParams
using ClimaLand
using SurfaceFluxes
import ClimaLand.Parameters as LP
using ClimaLand.InlandWater
using Dates
using ClimaLand.Domains: Point
using ClimaLand: prescribed_forcing_era5
using ClimaUtilities.TimeManager: ITime
using ClimaCore
FT = Float32
toml_dict = LP.create_toml_dict(FT)
default_params_from_toml = SlabLakeParameters(toml_dict)
longlat = FT.((50.6689, 41.9350))
z_sfc = FT(0)
domain = Point(; z_sfc, longlat)
surface_space = domain.space.surface;
start_date = DateTime(2008);
stop_date = start_date + Day(1)
forcing = prescribed_forcing_era5(
    start_date,
    stop_date,
    domain.space.surface,
    toml_dict,
    FT;
    use_lowres_forcing = true,
)
lake = InlandWater.SlabLakeModel(FT, domain, forcing, toml_dict)
@test ClimaLand.name(lake) == :lake
Y, p, cds = ClimaLand.initialize(lake)
@test propertynames(Y.lake) == prognostic_vars(lake)
@test propertynames(p.lake) == auxiliary_vars(lake)
@test InlandWater.get_turb_fluxes_type(FT, lake.boundary_conditions) ==
      NamedTuple{
    (:lhf, :shf, :vapor_flux, :∂lhf∂T, :∂shf∂T),
    Tuple{FT, FT, FT, FT, FT},
}
@test ClimaLand.surface_roughness_model(lake, Y, p) ==
      SurfaceFluxes.ConstantRoughnessParams{FT}(
    default_params_from_toml.z_0m,
    default_params_from_toml.z_0b,
)
@test ClimaLand.surface_height(lake, Y, p) == FT(0)
@test ClimaLand.surface_displacement_height(lake, Y, p) == FT(0)
@test all(parent(lake.inland_water_mask) .>= 0)
@test all(parent(lake.inland_water_mask) .<= 1)
t0 = ITime(0, Second(1), start_date)
ClimaLand.Simulations.set_lake_initial_conditions!(Y, p, t0, lake)
tmp = ClimaCore.Fields.zeros(surface_space)
evaluate!(tmp, forcing.atmos.T, t0)
expected_U =
    @. InlandWater.lake_energy_from_temperature(tmp, default_params_from_toml)
@test Y.lake.U == expected_U
set_initial_cache! = make_set_initial_cache(lake)
set_initial_cache!(p, Y, t0)
@test p.lake.T == tmp
@test p.lake.q_l == ClimaCore.Fields.zeros(surface_space)
@test all(parent(p.lake.albedo) .== default_params_from_toml.ice_albedo)
@test all(
    parent(
        (
            p.drivers.P_liq .+ p.drivers.P_snow .+
            p.lake.turbulent_fluxes.vapor_flux
        ) .+ p.lake.runoff,
    ) .≈ 0,
)
compute_exp_tendency! = make_compute_exp_tendency(lake)
dY = similar(Y)
compute_exp_tendency!(dY, Y, p, t0)
@test all(.~isnan.(parent(dY.lake.U)))
@test all(.~isinf.(parent(dY.lake.U)))
