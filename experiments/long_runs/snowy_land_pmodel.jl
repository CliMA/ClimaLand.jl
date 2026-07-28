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
# Timestep: 900 s
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
# If you want to run with prognostic (modeled) LAI from the optimal-LAI model
# (Zhou et al. 2025, `ZhouOptimalLAIModel`) instead of prescribed MODIS LAI,
# type `export PROGNOSTIC_LAI=""` in the terminal and run this script, or pass
# `PROGNOSTIC_LAI=""` as an environment variable on Buildkite. The default
# (unset) prescribes MODIS LAI.
const PROGNOSTIC_LAI = haskey(ENV, "PROGNOSTIC_LAI") ? true : false
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
lai_suffix = PROGNOSTIC_LAI ? "_opt_lai" : ""
# Optional RUN_TAG keeps concurrent config variants (e.g. z-cut vs default) in
# separate output directories instead of clashing on one root_path.
run_tag = get(ENV, "RUN_TAG", "")
tag_suffix = isempty(run_tag) ? "" : "_$(run_tag)"
# OUTPUT_ROOT relocates the (multi-GB) run output off the repo/home filesystem;
# default "." preserves the in-place behavior. On Derecho, point it at scratch.
output_root = get(ENV, "OUTPUT_ROOT", ".")
root_path = joinpath(
    output_root,
    "snowy_land_pmodel$(lai_suffix)$(tag_suffix)_longrun_$(device_suffix)",
)
diagnostics_outdir = joinpath(root_path, "global_diagnostics")
outdir =
    ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir)

function setup_model(
    ::Type{FT},
    start_date,
    stop_date,
    Δt,
    domain,
    toml_dict;
    prognostic_lai = false,
) where {FT}
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
    )
    forcing = (; atmos, radiation)

    prognostic_land_components = (:canopy, :lake, :snow, :soil, :soilco2)
    if prognostic_lai
        # The LandModel constructor uses the prognostic LAI model if no
        # prescribed LAI is passed.
        land = LandModel{FT}(
            forcing,
            toml_dict,
            domain,
            Δt;
            prognostic_land_components,
        )
    else
        # Prescribed LAI (default): read LAI from MODIS data.
        LAI = ClimaLand.Canopy.prescribed_lai_modis(
            surface_space,
            start_date,
            stop_date,
        )
        land = LandModel{FT}(
            forcing,
            LAI,
            toml_dict,
            domain,
            Δt;
            prognostic_land_components,
        )
    end
    return land
end

# If not LONGER_RUN, run for RUN_YEARS years (default 2); the forcing from 2008 is
# repeated. Set RUN_YEARS higher to allow a multi-year spin-up before scoring (e.g. 3 =
# 2-year spin-up + 1 scored year for the optimal-LAI running means). If LONGER run, run
# for 19 years with the correct forcing each year. Note that since the Northern
# hemisphere's winter season is defined as DJF, we simulate from and until the beginning
# of March so that a full season is included in seasonal metrics.
start_date = LONGER_RUN ? DateTime("2000-03-01") : DateTime("2008-03-01")
run_years = parse(Int, get(ENV, "RUN_YEARS", "2"))
stop_date = LONGER_RUN ? DateTime("2019-03-01") : start_date + Year(run_years)
Δt = 900.0
domain =
    ClimaLand.Domains.global_box_domain(FT; context, mask_threshold = FT(0.99))

# Optional ENV parameter overrides (mirrors experiments/integrated/era5/single_site.jl)
# so the global run can be launched with the optimal-LAI tuning knobs (online/band-pass
# f0, energy-cap z) without editing the calibrated toml. Unset ENV = calibrated default.
opt_overrides = Dict{String, FT}()
for (key, envname) in (
    ("optimal_lai_online_f0", "ONLINE_F0"),
    ("optimal_lai_f0_scale", "F0_SCALE"),
    ("optimal_lai_beta_in_A0", "BETA_IN_A0"),
    ("optimal_lai_online_vpd_gs", "ONLINE_VPD_GS"),
    ("optimal_lai_online_gsl", "ONLINE_GSL"),
    ("optimal_lai_online_c3c4", "ONLINE_C3C4"),
    ("optimal_lai_z", "OPT_Z"),
    ("optimal_lai_z_c4", "OPT_Z_C4"),
    ("optimal_lai_sigma", "OPT_SIGMA"),
    ("optimal_lai_sigma_c4", "OPT_SIGMA_C4"),
    ("optimal_lai_alpha", "OPT_ALPHA"),
    ("optimal_lai_z_a0", "Z_A0"),
)
    ev = get(ENV, envname, "")
    isempty(ev) || (opt_overrides[key] = parse(FT, ev))
end

override_files = String[]
UNCALIBRATED && push!(override_files, "toml/uncalibrated_parameters.toml")
if !isempty(opt_overrides)
    mkpath(root_path)
    override_path = joinpath(root_path, "override_params.toml")
    open(override_path, "w") do io
        for (k, v) in opt_overrides
            println(io, "[\"$k\"]")
            println(io, "value = $v")
            println(io, "type = \"float\"")
        end
    end
    push!(override_files, override_path)
end
toml_dict =
    isempty(override_files) ? LP.create_toml_dict(FT) :
    LP.create_toml_dict(FT; override_files)

model = setup_model(
    FT,
    start_date,
    stop_date,
    Δt,
    domain,
    toml_dict;
    prognostic_lai = PROGNOSTIC_LAI,
)
simulation = LandSimulation(start_date, stop_date, Δt, model; outdir)
@info "Run: Global Soil-Canopy-Snow Model"
@info "LAI: $(PROGNOSTIC_LAI ? "prognostic (ZhouOptimalLAIModel)" : "prescribed (MODIS)")"
@info "Resolution: $(domain.nelements)"
@info "Timestep: $Δt s"
@info "Start Date: $start_date"
@info "Stop Date: $stop_date"
CP.log_parameter_information(toml_dict, joinpath(root_path, "parameters.toml"))
ClimaLand.Simulations.solve!(simulation)

LandSimVis.make_annual_timeseries(simulation; savedir = root_path)
LandSimVis.make_heatmaps(simulation; savedir = root_path, date = stop_date)
LandSimVis.make_leaderboard_plots(
    simulation;
    savedir = root_path,
    leaderboard_data_sources = ["ERA5", "ILAMB", "FlagshipCarbonMetrics"],
)

if LONGER_RUN
    include("../ilamb/ilamb_conversion.jl")
    make_compatible_with_ILAMB(
        joinpath(root_path, "global_diagnostics", "output_active"),
        joinpath(root_path, "global_diagnostics", "ILAMB_diagnostics"),
    )
end
