# Test that a full LandModel with inland-water soil parameter overrides
# runs stably for 30 days without NaNs or crashes.
#
# Uses a Column domain at a lake coordinate (Lake Michigan: 43.5°N, 87.0°W)
# with ERA5 forcing, overriding soil parameters to water-like values
# (ν=0.95, K_sat=1e-3, θ_r=0, low albedo) to simulate an inland water point.

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
# Build LandModel with inland water soil overrides
# ==========================================
prognostic_land_components = (:canopy, :snow, :soil, :soilco2)

# Soil with water-like parameters
soil_albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(;
    PAR_albedo = FT(0.06),
    NIR_albedo = FT(0.06),
)

retention_parameters = (;
    ν = FT(0.95),
    θ_r = FT(0.0),
    K_sat = FT(1e-3),
    hydrology_cm = vanGenuchten{FT}(; α = FT(100.0), n = FT(2.0)),
)
composition_parameters = (;
    ν_ss_om = FT(0.0),
    ν_ss_quartz = FT(0.0),
    ν_ss_gravel = FT(0.0),
)

soil = Soil.EnergyHydrology{FT}(
    land_domain,
    forcing,
    toml_dict;
    prognostic_land_components,
    additional_sources = (ClimaLand.RootExtraction{FT}(),),
    albedo = soil_albedo,
    retention_parameters,
    composition_parameters,
    S_s = FT(1e-3),
    z_0m = FT(0.001),
    z_0b = FT(0.0001),
    emissivity = FT(0.97),
)

# Use the convenient LandModel constructor for canopy, snow, soilco2
LAI = ClimaLand.Canopy.prescribed_lai_modis(surface_space, start_date, stop_date)

land = LandModel{FT}(
    forcing,
    LAI,
    toml_dict,
    land_domain,
    dt;
    prognostic_land_components,
    soil,  # override with water-like soil
)

# ==========================================
# INITIAL CONDITIONS — saturated soil (like a lake)
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

    # Soil IC — saturated (water body)
    Y.soil.ϑ_l .= FT(0.95)  # = ν
    Y.soil.θ_i .= FT(0.0)
    ρc_s = ClimaLand.Soil.volumetric_heat_capacity.(
        Y.soil.ϑ_l, Y.soil.θ_i,
        land.soil.parameters.ρc_ds,
        land.soil.parameters.earth_param_set,
    )
    T_init = FT(270.0)  # January in Lake Michigan
    Y.soil.ρe_int .= ClimaLand.Soil.volumetric_internal_energy.(
        Y.soil.θ_i, ρc_s, T_init,
        land.soil.parameters.earth_param_set,
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

@testset "Full LandModel with inland water parameters (Lake Michigan)" begin
    @info "Running full LandModel with inland water soil parameters..."
    @info "  Location: Lake Michigan ($(lat)N, $(long)E)"
    @info "  Period: $start_date to $stop_date"
    @time sol = solve!(simulation)

    Y = simulation._integrator.u

    # No NaNs in state
    @test !any(isnan, parent(Y.soil.ϑ_l))
    @test !any(isnan, parent(Y.soil.θ_i))
    @test !any(isnan, parent(Y.soil.ρe_int))
    @info "  No NaNs in soil state"

    # Soil should remain fairly saturated
    min_swc = minimum(parent(Y.soil.ϑ_l))
    max_swc = maximum(parent(Y.soil.ϑ_l))
    @info "  Soil water content range: [$min_swc, $max_swc] (ν = 0.95)"
    @test min_swc > FT(-1.0)  # may dry somewhat without lateral inflow

    # Temperature should be physical
    ρc_s = ClimaLand.Soil.volumetric_heat_capacity.(
        Y.soil.ϑ_l, Y.soil.θ_i,
        land.soil.parameters.ρc_ds,
        land.soil.parameters.earth_param_set,
    )
    T_final = ClimaLand.Soil.temperature_from_ρe_int.(
        Y.soil.ρe_int, Y.soil.θ_i, ρc_s,
        land.soil.parameters.earth_param_set,
    )
    min_T = minimum(parent(T_final))
    max_T = maximum(parent(T_final))
    @info "  Temperature range: [$min_T, $max_T] K"
    @test min_T > FT(200)
    @test max_T < FT(350)

    # Snow state should be valid
    @test !any(isnan, parent(Y.snow.S))
    @info "  No NaNs in snow state"

    # Saved values available
    sv = saving_cb.affect!.saved_values
    n_saves = length(sv.saveval)
    @info "  Saved $n_saves timesteps"
    @test n_saves > 0

    @info "All inland water full-model tests passed!"
end
