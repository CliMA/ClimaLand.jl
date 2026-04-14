# # Global run of land model

# The code sets up and runs ClimaLand v1, which
# includes soil, canopy, and snow, on a spherical domain,
# using ERA5 data as forcing. In this simulation, we have
# turned lateral flow off because horizontal boundary conditions and the
# land/sea mask are not yet supported by ClimaCore.

# Simulation Setup
# Number of spatial elements: 101 in horizontal, 15 in vertical
# Soil depth: 50 m
# Simulation duration: 730 d
# Timestep: 450 s
# Timestepper: ARS111
# Fixed number of iterations: 3
# Jacobian update: every new Newton iteration
# Atmos forcing update: every 3 hours

import ClimaComms
ClimaComms.@import_required_backends
using ClimaUtilities.ClimaArtifacts
import ClimaUtilities.TimeManager: ITime, date

import ClimaDiagnostics
import ClimaUtilities

import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP
using ClimaCore
using ClimaLand
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!

using Dates

using CairoMakie, GeoMakie, ClimaAnalysis
import ClimaLand.LandSimVis as LandSimVis

const FT = Float64;
# If you want to do a very long run locally, you can enter `export
# LONGER_RUN=""` in the terminal and run this script. If you want to do a very
# long run on Buildkite manually, then make a new build and pass `LONGER_RUN=""`
# as an environment variable. In both cases, the value of `LONGER_RUN` does not
# matter.
const LONGER_RUN = haskey(ENV, "LONGER_RUN") ? true : false
# If you want to do run the simulation with uncalibrated parameters, type
# `export UNCALIBRATED=""` in the terminal and run this script, or
# pass `UNCALIBRATED=""` as an environment variable on buildkite.
const UNCALIBRATED = haskey(ENV, "UNCALIBRATED") ? true : false
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
root_path = "snowy_land_pmodel_longrun_$(device_suffix)"
diagnostics_outdir = joinpath(root_path, "global_diagnostics")
outdir =
    ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir)

function setup_model(
    ::Type{FT},
    start_date,
    stop_date,
    Δt,
    domain,
    toml_dict,
) where {FT}
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    surface_space = domain.space.surface
    # Forcing data - high resolution
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        start_date,
        stop_date,
        surface_space,
        toml_dict,
        FT;
        max_wind_speed = 25.0,
        context,
        use_lowres_forcing = true, #TODO: remove if run on CliMA
    )
    forcing = (; atmos, radiation)

    # Read in LAI from MODIS data
    LAI = ClimaLand.Canopy.prescribed_lai_modis(
        surface_space,
        start_date,
        stop_date,
    )

    ground = ClimaLand.PrognosticGroundConditions{FT}()
    canopy_forcing = (; atmos, radiation, ground)
    prognostic_land_components = (:canopy, :lake, :snow, :soil, :soilco2)

    # Construct the P model manually since it is not a default
    photosynthesis = PModel{FT}(domain, toml_dict)
    conductance = PModelConductance{FT}(toml_dict)
    # Use the soil moisture stress function based on soil moisture only
    soil_moisture_stress =
        ClimaLand.Canopy.PiecewiseMoistureStressModel{FT}(domain, toml_dict)
    biomass =
        ClimaLand.Canopy.PrescribedBiomassModel{FT}(domain, LAI, toml_dict)
    canopy = ClimaLand.Canopy.CanopyModel{FT}(
        surface_domain,
        canopy_forcing,
        LAI,
        toml_dict;
        prognostic_land_components,
        photosynthesis,
        conductance,
        soil_moisture_stress,
        biomass,
    )

    # Snow model setup
    # Set β = 0 in order to regain model without density dependence
    α_snow = Snow.ZenithAngleAlbedoModel(toml_dict)
    horz_degree_res =
        sum(ClimaLand.Domains.average_horizontal_resolution_degrees(domain)) / 2 # mean of resolution in latitude and longitude, in degrees
    #horz_degree_res = FT(1.0)  # placeholder for column run
    scf = Snow.WuWuSnowCoverFractionModel(toml_dict, horz_degree_res)
    snow = Snow.SnowModel(
        FT,
        surface_domain,
        forcing,
        toml_dict,
        Δt;
        prognostic_land_components,
        α_snow,
        scf,
    )

    # ── Soil model with custom albedo ───────────────────────────────────
    # Load silt/clay from custom SoilGrids file for OfflineLinearSoilAlbedo
    raw_comp = Soil.soil_composition_parameters(
        domain.space.subsurface,
        FT;
        load_silt_clay = true,
        path = "/resnick/groups/esm/xwu/ClimaArtifacts/soilgrids/soilgrids_lowres/soil_solid_vol_fractions_soilgrids_lowres.nc",
    )
    #ν_ss_silt_val = FT(parent(ClimaLand.Domains.top_center_to_surface(raw_comp.ν_ss_silt))[1])
    #ν_ss_clay_val = FT(parent(ClimaLand.Domains.top_center_to_surface(raw_comp.ν_ss_clay))[1])
    
    ν_ss_silt_field = raw_comp.ν_ss_silt   # subsurface Field
    ν_ss_clay_field = raw_comp.ν_ss_clay   # subsurface Field

    # WSA
    α_soil = Soil.OfflineLinearSoilAlbedo{FT}(;
        intercept_BSA_vis = FT(-2.917),
        intercept_BSA_nir = FT(-1.384),
        coef_ν_BSA_vis = FT(-0.142),
        coef_ν_BSA_nir = FT(-0.142),
        coef_om_BSA_vis = FT(-0.244),
        coef_om_BSA_nir = FT(-0.032),
        coef_cf_BSA_vis = FT(0.106),
        coef_cf_BSA_nir = FT(0.032),
        coef_sand_BSA_vis = FT(0.103),
        coef_sand_BSA_nir = FT(0.124),
        coef_clay_BSA_vis = FT(0.235),
        coef_clay_BSA_nir = FT(0.267),
        coef_silt_BSA_vis = FT(0.262),
        coef_silt_BSA_nir = FT(0.100),
        coef_n_BSA_vis = FT(3.272),
        coef_n_BSA_nir = FT(4.211),
        coef_θ_BSA_vis = FT(-0.026),
        coef_θ_BSA_nir = FT(-0.016),
        ν_ss_silt = ν_ss_silt_field,
        ν_ss_clay = ν_ss_clay_field,
        α_min = FT(0.00),
        α_max = FT(1.00),
    )    

    soil = Soil.EnergyHydrology{FT}(
        domain,
        forcing,
        toml_dict;
        prognostic_land_components,
        additional_sources = (ClimaLand.RootExtraction{FT}(),),
        albedo = α_soil,
    )

    # Construct the land model with all default components except for snow
    land = LandModel{FT}(
        forcing,
        LAI,
        toml_dict,
        domain,
        Δt;
        prognostic_land_components,
        snow,
        canopy,
        soil,
    )
    return land
end

# If not LONGER_RUN, run for 2 years; note that the forcing from 2008 is repeated.
# If LONGER run, run for 19 years, with the correct forcing each year.
# Note that since the Northern hemisphere's winter season is defined as DJF,
# we simulate from and until the beginning of
# March so that a full season is included in seasonal metrics.
start_date = LONGER_RUN ? DateTime("2000-03-01") : DateTime("2008-03-01")
stop_date = LONGER_RUN ? DateTime("2019-03-01") : DateTime("2010-03-01")
Δt = 450.0
domain =
    ClimaLand.Domains.global_box_domain(FT; context, mask_threshold = FT(0.99))
#longlat = FT.((25, 25))
#zlim = FT.((-15, 0))
#nelements = 15
#dz_tuple = FT.((3, 0.05))
#domain = ClimaLand.Domains.Column(; zlim, longlat, nelements, dz_tuple)

if UNCALIBRATED
    override_params_path = "toml/uncalibrated_parameters.toml"
    toml_dict = LP.create_toml_dict(FT, override_files = [override_params_path])
else
    toml_dict = LP.create_toml_dict(FT)
end

# Add albedo diagnostics before LandSimulation
saveat_albedo = Second(3600)  # hourly saves to avoid memory issues
saving_cb = ClimaLand.NonInterpSavingCallback(start_date, stop_date, saveat_albedo)
sv = saving_cb.affect!.saved_values

model = setup_model(FT, start_date, stop_date, Δt, domain, toml_dict)
simulation = LandSimulation(start_date, stop_date, Δt, model; outdir)
# Modify LandSimulation to include the callback
#simulation = LandSimulation(start_date, stop_date, Δt, model; outdir, user_callbacks = (saving_cb,))
@info "Run: Global Soil-Canopy-Snow Model"
@info "Resolution: $(domain.nelements)"
@info "Timestep: $Δt s"
@info "Start Date: $start_date"
@info "Stop Date: $stop_date"
CP.log_parameter_information(toml_dict, joinpath(root_path, "parameters.toml"))
ClimaLand.Simulations.solve!(simulation)

LandSimVis.make_annual_timeseries(simulation; savedir = root_path)
LandSimVis.make_heatmaps(simulation; savedir = root_path, date = stop_date)
LandSimVis.make_leaderboard_plots(simulation; savedir = root_path)

if LONGER_RUN
    include("../ilamb/ilamb_conversion.jl")
    make_compatible_with_ILAMB(
        joinpath(root_path, "global_diagnostics", "output_active"),
        joinpath(root_path, "global_diagnostics", "ILAMB_diagnostics"),
    )
end

# ── Column visualization ────────────────────────────────────────────────
#using CairoMakie

#dw = simulation.diagnostics[1].output_writer
#diag_dict = dw.dict
#diag_keys = collect(keys(diag_dict))
#@info "Available diagnostics: $diag_keys"

#function extract_timeseries(diag_dict, short_name)
#    if !haskey(diag_dict, short_name)
#        @warn "Diagnostic '$short_name' not found"
#        return nothing, nothing
#    end
#    od = diag_dict[short_name]
#    times = collect(keys(od))
#    days = [Float64(t.counter) / 86400 for t in times]
#    vals = [parent(od[t])[end] for t in times]
#    return days, vals
#end

#mkpath(root_path)

# Figure 1: Surface fluxes
#fig1 = Figure(size = (1600, 1200), fontsize = 20)
#for (i, (sn, yl)) in enumerate([
#    ("lhf_1M_average", "Latent Heat (W/m²)"),
#    ("lwu_1M_average", "LW Up (W/m²)"),
#    ("swd_1M_average", "SW Down (W/m²)"),
#    ("lwd_1M_average", "LW Down (W/m²)"),
#])
#    r, c = divrem(i - 1, 2) .+ (1, 1)
#    days, vals = extract_timeseries(diag_dict, sn)
#    if days !== nothing
#        ax = Axis(fig1[r, c], ylabel = yl, xlabel = "Days")
#        lines!(ax, days, vals)
#    end
#end
#Label(fig1[0, :], "Column ($(longlat[1])°E, $(longlat[2])°N): Surface Fluxes", fontsize = 24)
#CairoMakie.save(joinpath(root_path, "column_surface_fluxes.png"), fig1)

# Figure 2: State variables
#fig2 = Figure(size = (1600, 1200), fontsize = 20)
#for (i, (sn, yl)) in enumerate([
#    ("swc_1M_average", "Soil Water Content"),
#    ("swe_1M_average", "SWE (m)"),
#    ("snd_1M_average", "Snow Depth (m)"),
#    ("tr_1M_average", "Transpiration"),
#])
#    r, c = divrem(i - 1, 2) .+ (1, 1)
#    days, vals = extract_timeseries(diag_dict, sn)
#    if days !== nothing
#        ax = Axis(fig2[r, c], ylabel = yl, xlabel = "Days")
#        lines!(ax, days, vals)
#    end
#end
#Label(fig2[0, :], "Column ($(longlat[1])°E, $(longlat[2])°N): State Variables", fontsize = 24)
#CairoMakie.save(joinpath(root_path, "column_state_variables.png"), fig2)

#@info "Saved plots to $root_path"

# ── Figure: Soil Albedo ─────────────────────────────────────────────────
#sv_times = Dates.value.(Second.(sv.t .- sv.t[1]))
#fig_alb = Figure(size = (1200, 600), fontsize = 20)
#ax = Axis(fig_alb[1, 1], ylabel = "Albedo", xlabel = "Days",
#    title = "Column ($(longlat[1])°E, $(longlat[2])°N): Soil Albedo")
#lines!(ax, sv_times ./ 86400,
#    [parent(sv.saveval[k].soil.PAR_albedo)[1] for k in 1:length(sv_times)],
#    label = "PAR (VIS)")
#lines!(ax, sv_times ./ 86400,
#    [parent(sv.saveval[k].soil.NIR_albedo)[1] for k in 1:length(sv_times)],
#    label = "NIR")
#axislegend(ax, position = :rt)
#CairoMakie.save(joinpath(root_path, "column_soil_albedo.png"), fig_alb)
#@info "Saved column_soil_albedo.png"
