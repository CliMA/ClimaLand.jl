# Test that a full LandModel with a slab lake runs stably for 30 days
# without NaNs and produces physical lake temperatures.
#
# Uses a Column domain at a lake coordinate (Lake Michigan: 43.5°N, 87.0°W)
# with ERA5 forcing. The SlabLakeModel is a separate component in the LandModel,
# with fraction-based blending of soil boundary conditions.

import SciMLBase
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
import ClimaParams as CP
using Dates
using Test

using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaLand.InlandWater
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!

const FT = Float64
toml_dict = LP.create_toml_dict(FT)

# Lake Michigan coordinate
lat = FT(43.5)
long = FT(-87.0)

# Domain — same vertical structure as global runs
depth = FT(15.0)
dz_tuple = FT.((3.0, 0.05))
land_domain = Column(;
    zlim = (-depth, FT(0)),
    nelements = 15,
    dz_tuple,
    longlat = (long, lat),
)
surface_space = land_domain.space.surface

# Simulation period — 30 days in 2008 (matches lowres ERA5 data)
dt = FT(450)
start_date = DateTime(2008)
stop_date = start_date + Day(30)

# ERA5 forcing
atmos, radiation = ClimaLand.prescribed_forcing_era5(
    start_date,
    stop_date,
    surface_space,
    toml_dict,
    FT;
    use_lowres_forcing = true,
)
forcing = (; atmos, radiation)

# ==========================================
# Build LandModel with slab lake
# ==========================================
prognostic_land_components = (:canopy, :snow, :soil, :soilco2)

# Inland water mask: the whole column is inland water (value = 1)
inland_water_mask = ClimaCore.Fields.zeros(FT, surface_space) .+ FT(1.0)

# Use the convenient LandModel constructor which handles lake creation
LAI =
    ClimaLand.Canopy.prescribed_lai_modis(surface_space, start_date, stop_date)

land = LandModel{FT}(
    forcing,
    LAI,
    toml_dict,
    land_domain,
    dt;
    prognostic_land_components,
    inland_water_mask,
)

# ==========================================
# INITIAL CONDITIONS
# ==========================================
function set_ic!(Y, p, t0, land)
    atmos_model = land.soil.boundary_conditions.top.atmos
    if atmos_model isa ClimaLand.PrescribedAtmosphere
        ClimaLand.Simulations.evaluate!(p.drivers.T, atmos_model.T, t0)
    end
    # SoilCO2 IC
    Y.soilco2.CO2 .= FT(0.000412)
    Y.soilco2.O2_f .= FT(0.21)
    Y.soilco2.SOC .= FT(5.0)

    # Soil sediment IC — saturated, at lake temperature
    T_init = FT(277.0)  # 4°C — typical late-fall lake temperature
    Y.soil.ϑ_l .= land.soil.parameters.ν
    Y.soil.θ_i .= FT(0.0)
    ρc_s =
        ClimaLand.Soil.volumetric_heat_capacity.(
            Y.soil.ϑ_l,
            Y.soil.θ_i,
            land.soil.parameters.ρc_ds,
            land.soil.parameters.earth_param_set,
        )
    Y.soil.ρe_int .=
        ClimaLand.Soil.volumetric_internal_energy.(
            Y.soil.θ_i,
            ρc_s,
            T_init,
            land.soil.parameters.earth_param_set,
        )

    # Slab lake IC — initialize U from temperature
    lake_params = land.lake.parameters
    earth_param_set = lake_params.earth_param_set
    Y.lake.U .= InlandWater.lake_energy_from_temperature(
        T_init,
        lake_params,
        earth_param_set,
    )

    # Snow IC — none
    Y.snow.S .= FT(0)
    Y.snow.S_l .= FT(0)
    Y.snow.U .= FT(0)

    # Canopy IC
    if land.canopy.energy isa ClimaLand.Canopy.BigLeafEnergyModel
        Y.canopy.energy.T .= p.drivers.T
    end
    Y.canopy.hydraulics.ϑ_l.:1 .= land.canopy.hydraulics.parameters.ν
end

saving_cb = ClimaLand.NonInterpSavingCallback(start_date, stop_date, Second(dt))

simulation = LandSimulation(
    start_date,
    stop_date,
    dt,
    land;
    set_ic!,
    user_callbacks = (saving_cb,),
    diagnostics = nothing,
)

@testset "Slab lake LandModel at Lake Michigan (30-day)" begin
    @info "Running slab lake LandModel..."
    @info "  Location: Lake Michigan ($(lat)N, $(long)E)"
    @info "  Period: $start_date to $stop_date"
    @time sol = solve!(simulation)

    Y = simulation._integrator.u
    p = simulation._integrator.p

    # --- Lake component exists ---
    @test hasproperty(Y, :lake)
    @test hasproperty(Y.lake, :U)
    @test !any(isnan, parent(Y.lake.U))

    # --- No NaNs in soil sediment state ---
    @test !any(isnan, parent(Y.soil.ϑ_l))
    @test !any(isnan, parent(Y.soil.θ_i))
    @test !any(isnan, parent(Y.soil.ρe_int))

    # --- Lake temperature is physical ---
    T_lake = parent(p.lake.T)[1]
    @info "  T_lake = $T_lake K"
    @test T_lake > FT(250)   # not implausibly cold
    @test T_lake < FT(350)   # not implausibly hot

    # --- Lake liquid fraction is in [0, 1] ---
    q_l = parent(p.lake.q_l)[1]
    @info "  q_l = $q_l"
    @test FT(0) <= q_l <= FT(1)

    # --- Sediment temperatures are physical ---
    T_soil = parent(p.soil.T)
    min_T = minimum(T_soil)
    max_T = maximum(T_soil)
    @info "  Sediment temperature range: [$min_T, $max_T] K"
    @test min_T > FT(250)
    @test max_T < FT(350)

    # --- Surface energy fluxes are physical for a lake ---
    lhf = parent(p.lake.turbulent_fluxes.lhf)[1]
    shf = parent(p.lake.turbulent_fluxes.shf)[1]
    R_n = parent(p.lake.R_n)[1]
    sfc_energy = parent(p.lake.surface_energy_flux)[1]
    sed_heat = parent(p.lake.sediment_heat_flux)[1]
    @info "  Lake surface energy fluxes (W/m²):"
    @info "    LHF  = $lhf"
    @info "    SHF  = $shf"
    @info "    R_n  = $R_n"
    @info "    Lake surface energy flux = $sfc_energy"
    @info "    Lake–sediment heat flux  = $sed_heat"
    # Fluxes should be finite and within reasonable bounds for a winter lake
    @test isfinite(lhf)
    @test isfinite(shf)
    @test isfinite(R_n)
    @test abs(lhf) < FT(500)   # latent heat < 500 W/m²
    @test abs(shf) < FT(500)   # sensible heat < 500 W/m²
    @test abs(R_n) < FT(1000)  # net radiation < 1000 W/m²

    # --- Snow state is valid ---
    @test !any(isnan, parent(Y.snow.S))

    # --- Saved values available ---
    sv = saving_cb.affect!.saved_values
    n_saves = length(sv.saveval)
    @info "  Saved $n_saves timesteps"
    @test n_saves > 0
end
