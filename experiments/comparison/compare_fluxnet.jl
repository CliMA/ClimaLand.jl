# Compare the single-column `Column` build against a single-point `ColumnEnsemble`
# (multi-column `MultiColumnFiniteDifferenceSpace`) for a FULL integrated fluxnet
# simulation — the most feature-rich model in the repo. The model setup is copied
# from `run_fluxnet.jl`: an integrated `LandModel` of
#   - soil: `EnergyHydrology` with root extraction + surface runoff,
#   - soil CO2 biogeochemistry (`SoilCO2Model`, prognostic microbes),
#   - canopy: two-stream radiative transfer, Medlyn conductance, Farquhar
#     photosynthesis, plant hydraulics, big-leaf energy, prescribed biomass,
#   - snow (`SnowModel`),
# driven by real US-MOz fluxnet forcing and MODIS LAI.
#
# Both builds share identical vertical mesh, parameters, forcing, and initial
# condition, so column 1 of the ensemble must reproduce the single column. Every
# prognostic and cache field is compared (via `report_state_diffs` from
# field_rmse.jl) at t = 0 before any steps, then again after a single step.
#
# Run: julia --project=.buildkite experiments/comparison/compare_fluxnet.jl [SITE_ID]
#      (SITE_ID defaults to US-MOz)

import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
import ClimaParams as CP
using Dates

using ClimaLand
using ClimaLand.Domains: Column, ColumnEnsemble
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!, step!
import ClimaLand.FluxnetSimulations as FluxnetSimulations

include(joinpath(@__DIR__, "field_rmse.jl"))

const FT = Float64
toml_dict = LP.create_toml_dict(FT)
prognostic_land_components = (:canopy, :snow, :soil, :soilco2)

site_ID = isempty(ARGS) ? "US-MOz" : ARGS[1]
site_ID_val = FluxnetSimulations.replace_hyphen(site_ID)

# Site defaults: vertical mesh, location, tower height, and all soil/plant params.
(; dz_tuple, nelements, zmin, zmax) =
    FluxnetSimulations.get_domain_info(FT, Val(site_ID_val))
(; time_offset, lat, long) =
    FluxnetSimulations.get_location(FT, Val(site_ID_val))
(; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))
(;
    soil_ν,
    soil_K_sat,
    soil_S_s,
    soil_vg_n,
    soil_vg_α,
    θ_r,
    ν_ss_quartz,
    ν_ss_om,
    ν_ss_gravel,
    z_0m_soil,
    z_0b_soil,
    soil_ϵ,
    soil_α_PAR,
    soil_α_NIR,
    Ω,
    χl,
    α_PAR_leaf,
    λ_γ_PAR,
    τ_PAR_leaf,
    α_NIR_leaf,
    τ_NIR_leaf,
    ϵ_canopy,
    ac_canopy,
    g1,
    Drel,
    g0,
    Vcmax25,
    SAI,
    f_root_to_shoot,
    K_sat_plant,
    ψ63,
    Weibull_param,
    a,
    conductivity_model,
    retention_model,
    plant_ν,
    plant_S_s,
    rooting_depth,
    h_canopy,
) = FluxnetSimulations.get_parameters(FT, Val(site_ID_val))

dt = Float64(450) # 7.5 minutes

# Shared forcing (atmos + radiation) from the flux tower data; identical for both builds.
(start_date, stop_date) =
    FluxnetSimulations.get_data_dates(site_ID, time_offset)
(; atmos, radiation) = FluxnetSimulations.prescribed_forcing_fluxnet(
    site_ID,
    lat,
    long,
    time_offset,
    atmos_h,
    start_date,
    toml_dict,
    FT,
)
forcing = (; atmos, radiation)

# Maximum LAI is location-based and shared; computed once.
maxLAI = FluxnetSimulations.get_maxLAI_at_site(start_date, lat, long)
RAI = maxLAI * f_root_to_shoot

# Build the full integrated LandModel + simulation on the given domain. This is
# the run_fluxnet.jl model setup, parameterized by `land_domain` so we can build
# it identically on a Column and a single-point ColumnEnsemble.
function build_fluxnet_sim(land_domain)
    surface_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)
    soil_domain = land_domain

    soil_albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(;
        PAR_albedo = soil_α_PAR,
        NIR_albedo = soil_α_NIR,
    )
    retention_parameters = (;
        ν = soil_ν,
        θ_r,
        K_sat = soil_K_sat,
        hydrology_cm = vanGenuchten{FT}(; α = soil_vg_α, n = soil_vg_n),
    )
    composition_parameters = (; ν_ss_om, ν_ss_quartz, ν_ss_gravel)
    soil = Soil.EnergyHydrology{FT}(
        soil_domain,
        forcing,
        toml_dict;
        prognostic_land_components,
        additional_sources = (ClimaLand.RootExtraction{FT}(),),
        albedo = soil_albedo,
        runoff = ClimaLand.Soil.Runoff.SurfaceRunoff(),
        retention_parameters,
        composition_parameters,
        S_s = soil_S_s,
        z_0m = z_0m_soil,
        z_0b = z_0b_soil,
        emissivity = soil_ϵ,
    )

    co2_prognostic_soil = Soil.Biogeochemistry.PrognosticMet(soil.parameters)
    co2_drivers = Soil.Biogeochemistry.SoilDrivers(co2_prognostic_soil, atmos)
    soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(
        soil_domain,
        co2_drivers,
        toml_dict,
    )

    radiation_parameters = (;
        Ω,
        G_Function = CLMGFunction(χl),
        α_PAR_leaf,
        τ_PAR_leaf,
        α_NIR_leaf,
        τ_NIR_leaf,
    )
    radiative_transfer = Canopy.TwoStreamModel{FT}(
        surface_domain,
        toml_dict;
        radiation_parameters,
        ϵ_canopy,
    )
    conductance =
        Canopy.MedlynConductanceModel{FT}(surface_domain, toml_dict; g1)
    photosynthesis_parameters = (; fractional_c3 = FT(1), Vcmax25)
    photosynthesis =
        FarquharModel{FT}(surface_domain, toml_dict; photosynthesis_parameters)

    surface_space = land_domain.space.surface
    LAI = ClimaLand.Canopy.prescribed_lai_modis(
        surface_space,
        start_date,
        stop_date,
    )
    hydraulics = Canopy.PlantHydraulicsModel{FT}(
        surface_domain,
        toml_dict;
        ν = plant_ν,
        S_s = plant_S_s,
        conductivity_model,
        retention_model,
    )
    biomass = Canopy.PrescribedBiomassModel{FT}(;
        LAI,
        SAI,
        RAI,
        rooting_depth,
        height = h_canopy,
    )
    energy = Canopy.BigLeafEnergyModel{FT}(toml_dict; ac_canopy)
    ground = ClimaLand.PrognosticGroundConditions{FT}()
    canopy_forcing = (; atmos, radiation, ground)
    canopy = Canopy.CanopyModel{FT}(
        surface_domain,
        canopy_forcing,
        LAI,
        toml_dict;
        prognostic_land_components,
        radiative_transfer,
        photosynthesis,
        conductance,
        hydraulics,
        energy,
        biomass,
    )

    snow = Snow.SnowModel(
        FT,
        surface_domain,
        forcing,
        toml_dict,
        dt;
        prognostic_land_components,
    )

    lake = nothing
    land = LandModel{FT}(canopy, snow, soil, soilco2, lake)
    set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
        site_ID,
        start_date,
        time_offset,
        land,
    )
    return LandSimulation(
        start_date,
        stop_date,
        dt,
        land;
        set_ic!,
        updateat = dt,
        diagnostics = nothing,
    )
end

column_domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = dz_tuple,
    longlat = (long, lat),
)
ensemble_domain = ColumnEnsemble(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = dz_tuple,
    longlat = (long, lat),
)

sim_column = build_fluxnet_sim(column_domain)
sim_ensemble = build_fluxnet_sim(ensemble_domain)

# Per field, with a1 = flattened single-column field, a2 = flattened ensemble
# column, and n = length(a1):
#   rmse  = sqrt(sum((a1 .- a2) .^ 2) / n)   RMS difference (field units)
#   scale = sqrt(sum(a1 .^ 2) / n)           RMS magnitude of the reference a1
#   rel   = rmse / scale                     relative error (dimensionless;
#                                            0 if rmse = scale = 0, Inf if only scale = 0)
report_state_diffs(
    sim_column,
    sim_ensemble;
    col = 1,
    label = "t = 0 (no steps)",
)

step!(sim_column)
step!(sim_ensemble)

@assert sum(isnan.(sim_ensemble._integrator.u)) == 0

report_state_diffs(sim_column, sim_ensemble; col = 1, label = "after 1 step")

@info "Full integrated fluxnet LandModel ran on the multi-column \
       MultiColumnFiniteDifferenceSpace and reproduced the single Column."
