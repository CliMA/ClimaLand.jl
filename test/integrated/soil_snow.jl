using Test
import ClimaComms
ClimaComms.@import_required_backends
using ClimaLand
import ClimaLand
import ClimaParams as CP
import ClimaLand.Parameters as LP
using Dates
FT = Float32
toml_dict = LP.create_toml_dict(FT)
earth_param_set = LP.LandParameters(toml_dict)
domain =
    ClimaLand.Domains.global_domain(FT; nelements = (5, 15), apply_mask = false)
(atmos, radiation) = prescribed_analytic_forcing(FT; toml_dict)
forcing = (; atmos, radiation)
Δt = FT(180.0)
land_model = ClimaLand.SoilSnowModel{FT}(forcing, toml_dict, domain, Δt)

Y, p, coords = ClimaLand.initialize(land_model)
p_soil_alone = deepcopy(p)
for lsm_aux_var in (
    :excess_water_flux,
    :excess_heat_flux,
    :ground_heat_flux,
    :effective_soil_sfc_T,
    :sfc_scratch,
    :subsfc_scratch,
    :effective_soil_sfc_depth,
)
    @test lsm_aux_var ∈ propertynames(p)
end
@test typeof(land_model.soil.boundary_conditions.top) <:
      ClimaLand.Soil.AtmosDrivenFluxBC
@test land_model.soil.boundary_conditions.top.prognostic_land_components ==
      (:snow, :soil)
src = ClimaLand.SoilSublimationwithSnow{FT}()
@test src ∈ land_model.soil.sources
@test ClimaLand.get_drivers(land_model) == (atmos, radiation)
# Set initial conditions for a case with *no snow on ground*
function init_soil!(Y, z, params)
    ν = params.ν
    FT = eltype(ν)
    @. Y.soil.ϑ_l = ν / 2
    Y.soil.θ_i .= 0.05
    T = FT(273)
    ρc_s = @. ClimaLand.Soil.volumetric_heat_capacity(
        ν / 2,
        FT(0),
        params.ρc_ds,
        params.earth_param_set,
    )
    Y.soil.ρe_int =
        ClimaLand.Soil.volumetric_internal_energy.(
            FT(0),
            ρc_s,
            T,
            params.earth_param_set,
        )
end
function init_snow!(Y, S, params)
    Y.snow.S .= S
    Y.snow.S_l .= 0
    @. Y.snow.U =
        ClimaLand.Snow.energy_from_q_l_and_swe(Y.snow.S, 0.0f0, params)
end

t = Float64(0)
init_soil!(Y, coords.subsurface.z, land_model.soil.parameters)
init_snow!(Y, 0.0f0, land_model.snow.parameters)

set_initial_cache! = make_set_initial_cache(land_model)
set_initial_cache!(p, Y, t)
# Make sure snow boundary fluxes are zero
@test all(parent(p.snow.total_energy_flux) .≈ 0)
@test all(parent(p.snow.total_water_flux) .≈ 0)
# Make sure the boundary conditions match bare soil result
set_soil_initial_cache! = make_set_initial_cache(land_model.soil)
set_soil_initial_cache!(p_soil_alone, Y, t)
@test p.soil.top_bc == p_soil_alone.soil.top_bc
dY_soil_snow = deepcopy(Y) .* 0
dY_soil_alone = deepcopy(Y) .* 0
ClimaLand.source!(
    dY_soil_alone,
    ClimaLand.Soil.SoilSublimation{FT}(),
    Y,
    p_soil_alone,
    land_model.soil,
)
ClimaLand.source!(dY_soil_snow, src, Y, p, land_model.soil)
@test dY_soil_alone.soil == dY_soil_snow.soil


# Repeat now with snow cover fraction of 1
init_snow!(Y, 1.0f0, land_model.snow.parameters)
set_initial_cache!(p, Y, t)
# Make sure the boundary conditions for soil are correct since that is what the LandHydrology methods
# affect
@test all(
    parent(p.soil.top_bc.water) .≈
    parent(p.excess_water_flux .+ p.snow.water_runoff),
)
@test all(
    parent(p.soil.top_bc.heat) .≈ parent(
        p.excess_heat_flux .+ p.snow.snow_cover_fraction .* p.ground_heat_flux,
    ),
)
dY_soil_snow = deepcopy(Y) .* 0
ClimaLand.source!(dY_soil_snow, src, Y, p, land_model.soil)
@test all(parent(dY_soil_snow.soil.θ_i) .≈ 0)

# Make sure soil boundary flux method also worked
G = deepcopy(p.ground_heat_flux)
p_snow_alone = deepcopy(p)
ClimaLand.Snow.snow_boundary_fluxes!(
    land_model.snow.boundary_conditions,
    Val((:snow,)),
    land_model.snow,
    Y,
    p_snow_alone,
    t,
)
ClimaLand.Snow.snow_boundary_fluxes!(
    land_model.snow.boundary_conditions,
    Val((:snow, :soil)),
    land_model.snow,
    Y,
    p,
    t,
)
@test all(
    parent(p_snow_alone.snow.total_energy_flux .- p.snow.total_energy_flux) .≈
    parent(p.snow.snow_cover_fraction .* p.ground_heat_flux),
)
