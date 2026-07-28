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
# OUTPUT_ROOT relocates the (multi-GB) output off the repo/home filesystem;
# default "." keeps the in-place behavior. On an HPC, point it at scratch.
output_root = get(ENV, "OUTPUT_ROOT", ".")
root_path = joinpath(
    output_root,
    "snowy_land_pmodel$(lai_suffix)_longrun_$(device_suffix)",
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

# If not LONGER_RUN, run for 2 years; note that the forcing from 2008 is repeated.
# If LONGER run, run for 19 years, with the correct forcing each year.
# Note that since the Northern hemisphere's winter season is defined as DJF,
# we simulate from and until the beginning of
# March so that a full season is included in seasonal metrics.
start_date = LONGER_RUN ? DateTime("2000-03-01") : DateTime("2008-03-01")
# RUN_YEARS (default 2) sets the non-LONGER_RUN length. Raise it to spin the
# optimal-LAI trailing totals up to equilibrium: they have a `tau_long_term`
# memory (2 years by default), so a 10-year run is roughly what it takes to
# produce an initial condition that no longer reflects the starting artifact.
run_years = parse(Int, get(ENV, "RUN_YEARS", "2"))
stop_date =
    LONGER_RUN ? DateTime("2019-03-01") : start_date + Dates.Year(run_years)
Δt = 900.0
domain =
    ClimaLand.Domains.global_box_domain(FT; context, mask_threshold = FT(0.99))

if UNCALIBRATED
    override_params_path = "toml/uncalibrated_parameters.toml"
    toml_dict = LP.create_toml_dict(FT, override_files = [override_params_path])
else
    toml_dict = LP.create_toml_dict(FT)
end

model = setup_model(
    FT,
    start_date,
    stop_date,
    Δt,
    domain,
    toml_dict;
    prognostic_lai = PROGNOSTIC_LAI,
)
# With prognostic LAI, also write the six fields an optimal-LAI initial
# condition is built from (see `set_canopy_component_initial_conditions!`), so a
# spun-up run can be turned back into the `optimal_lai_inputs` artifact.
diagnostics = if PROGNOSTIC_LAI
    ClimaLand.Diagnostics.default_diagnostics(
        model,
        start_date,
        outdir;
        output_vars = union(
            ClimaLand.Diagnostics.get_short_diagnostics(model),
            ["lai", "a0a", "pra", "olf0", "olvpd", "olgsl"],
        ),
    )
else
    ClimaLand.Diagnostics.default_diagnostics(model, start_date, outdir)
end
simulation =
    LandSimulation(start_date, stop_date, Δt, model; outdir, diagnostics)
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
