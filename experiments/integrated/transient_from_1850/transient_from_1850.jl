# TRENDY-style transient runs of the global LandModel with prognostic LAI and carbon pools,
# from a spun-up 1850 state (see README.md).
#
# A run covers the nominal years START_YEAR to END_YEAR. Nominal year n is forced by the
# forcing year f(n): the recycled 1901-1920 cycle, f(n) = 1901 + mod(n - 1850, 20), before
# 1901 or when the scenario recycles the climate, and f(n) = n otherwise. Consecutive
# nominal years with consecutive forcing years form a segment, run as one LandSimulation
# from the state at the end of the previous year.
#
# At the end of each year the state is saved to <RUN_DIR>/states/state_<n>.hdf5 and a row
# of global carbon fluxes and stocks (Pg C, land north of 60°S) is appended to
# <RUN_DIR>/annual_budget.csv. A rerun resumes after the last saved year, and a job stops
# at a year boundary before WALLTIME_HOURS, so long runs can be chained.
#
# Options (environment variables), with defaults:
#   SCENARIO=S2          S2: CO2 and climate; S1: CO2, recycled climate; S0: 1850 CO2,
#                        recycled climate (the control)
#   START_YEAR=1850, END_YEAR=2023
#   ARTIFACT_DIR         the folder of the climaland_1850_spunup_ic artifact (the 1850 state
#                        and the CO2 series), until it is registered in Artifacts.toml
#   IC=<ARTIFACT_DIR>/climaland_1850_spunup_ic.nc, or a checkpoint (.hdf5) of an earlier
#                        run to branch from
#   CO2_FILE=<ARTIFACT_DIR>/co2_annual_1700_2025.csv
#   FORCING=crujra       or era5 (from 1979)
#   RUN_DIR=transient_from_1850_<SCENARIO> (in the working directory)
#   WALLTIME_HOURS=12, MAX_YEARS (no limit), OUTPUT_VARS (comma-separated; default below)
#
#   julia --project=.buildkite experiments/integrated/transient_from_1850/transient_from_1850.jl
import ClimaComms
ClimaComms.@import_required_backends
import ClimaCore
import ClimaUtilities
import ClimaUtilities.TimeManager: date
import ClimaParams as CP
using ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation
using Dates

include(joinpath(@__DIR__, "co2_series.jl"))
include(joinpath(@__DIR__, "carbon_budget.jl"))
include(joinpath(@__DIR__, "initial_state.jl"))

const FT = Float64 # a year's change of SOC or wood is below Float32 resolution
const Δt = 900.0
const DEFAULT_OUTPUT_VARS = [
    "gpp",
    "nee",
    "er",
    "hr",
    "ra",
    "lai",
    "cveg",
    "soc_int",
    "ftr",
    "lhf",
    "shf",
    "swu",
    "lwu",
    "et",
    "tair",
    "precip",
    "swd",
]

env(name, default) = get(ENV, name, default)
const SCENARIO = env("SCENARIO", "S2")
SCENARIO in ("S0", "S1", "S2") || error("unknown SCENARIO $SCENARIO")
const START_YEAR = parse(Int, env("START_YEAR", "1850"))
const END_YEAR = parse(Int, env("END_YEAR", "2023"))
const ARTIFACT_DIR = env("ARTIFACT_DIR", "")
const IC = env("IC", joinpath(ARTIFACT_DIR, "climaland_1850_spunup_ic.nc"))
const FORCING = env("FORCING", "crujra")
FORCING in ("crujra", "era5") || error("unknown FORCING $FORCING")
const CO2_FILE =
    env("CO2_FILE", joinpath(ARTIFACT_DIR, "co2_annual_1700_2025.csv"))
const RUN_DIR = abspath(env("RUN_DIR", "transient_from_1850_$(SCENARIO)"))
const WALLTIME_HOURS = parse(Float64, env("WALLTIME_HOURS", "12"))
const MAX_YEARS =
    haskey(ENV, "MAX_YEARS") ? parse(Int, ENV["MAX_YEARS"]) : typemax(Int)
const OUTPUT_VARS =
    haskey(ENV, "OUTPUT_VARS") ? String.(split(ENV["OUTPUT_VARS"], ',')) :
    DEFAULT_OUTPUT_VARS
const STATE_DIR = joinpath(RUN_DIR, "states")
const BUDGET_CSV = joinpath(RUN_DIR, "annual_budget.csv")
const JOB_START = time()

const RECYCLE = SCENARIO in ("S0", "S1")
const CO2_FIXED_YEAR = SCENARIO == "S0" ? 1850 : nothing
forcing_year(n) = (RECYCLE || n < 1901) ? 1901 + mod(n - 1850, 20) : n
# CRU-JRA ends at 2023-12-31 18:00 (6-hourly), with no 2024 file: a run to 2023 stops
# there, and its last year has 364.75 days
const LAST_FORCING_TIME =
    FORCING == "crujra" ? DateTime(2023, 12, 31, 18) : DateTime(2100)

"""Consecutive nominal years from `n0` to `n1` with consecutive forcing years."""
function segment_end(n0, n1)
    n = n0
    while n < n1 && forcing_year(n + 1) == forcing_year(n) + 1
        n += 1
    end
    return n
end

state_path(n) = joinpath(STATE_DIR, "state_$(n).hdf5")

function last_saved_year()
    isdir(STATE_DIR) || return nothing
    years = Int[]
    for f in readdir(STATE_DIR)
        m = match(r"^state_(\d+)\.hdf5$", f)
        isnothing(m) || push!(years, parse(Int, m[1]))
    end
    return isempty(years) ? nothing : maximum(years)
end

function save_state(simulation, n)
    tmp = mktempdir(STATE_DIR)
    ClimaLand.save_checkpoint(
        simulation._integrator.u,
        simulation._integrator.t,
        tmp;
        model = simulation.model,
    )
    mv(only(readdir(tmp; join = true)), state_path(n); force = true)
    rm(tmp; recursive = true)
end

function make_land(start_date, stop_date, domain, toml_dict, co2)
    surface_space = domain.space.surface
    context = ClimaComms.context(surface_space)
    c_co2 = co2_time_varying_input(co2)
    forcing = if FORCING == "crujra"
        ClimaLand.prescribed_forcing_crujra(
            start_date,
            stop_date,
            surface_space,
            toml_dict,
            FT;
            c_co2,
            context,
            fill_missing = true,
            zenith_scaled_shortwave = true,
        )
    else
        ClimaLand.prescribed_forcing_era5(
            start_date,
            stop_date,
            surface_space,
            toml_dict,
            FT;
            max_wind_speed = 25.0,
            c_co2,
            context,
        )
    end
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    biomass = ClimaLand.Canopy.PrognosticCarbonModel{FT}(
        ClimaLand.Canopy.ZhouOptimalLAIModel{FT}(surface_domain, toml_dict),
        toml_dict,
    )
    return LandModel{FT}(
        forcing,
        toml_dict,
        domain,
        Δt;
        prognostic_land_components = (:canopy, :lake, :snow, :soil, :soilco2),
        biomass,
    )
end

"""Writes the soil CO2 and carbon pool parameters as built, once per run directory."""
function log_parameters(land, toml_dict)
    path = joinpath(RUN_DIR, "parameters_used.txt")
    isfile(path) && return nothing
    lines = String[]
    for (label, params) in (
        ("soilco2", land.soilco2.parameters),
        ("carbon pools", land.canopy.biomass.parameters),
    )
        for name in fieldnames(typeof(params))
            value = getfield(params, name)
            value isa Number && push!(lines, "$label.$name = $value")
        end
    end
    write(path, join(lines, "\n") * "\n")
    CP.log_parameter_information(
        toml_dict,
        joinpath(RUN_DIR, "parameters.toml"),
    )
end

"""Runs the nominal years n0..n1 of one segment; returns the last year completed."""
function run_segment(n0, n1, domain, toml_dict, ic, seconds_per_year)
    f0 = forcing_year(n0)
    start_date = DateTime(f0)
    stop_date = min(DateTime(f0 + n1 - n0 + 1), LAST_FORCING_TIME)
    co2 =
        isnothing(CO2_FIXED_YEAR) ?
        AnnualCO2Series(CO2_FILE; year_offset = n0 - f0) :
        AnnualCO2Series(CO2_FILE; fixed_year = CO2_FIXED_YEAR)
    land = make_land(start_date, stop_date, domain, toml_dict, co2)
    outdir = joinpath(RUN_DIR, "output", "seg_$(n0)_f$(f0)")
    mkpath(outdir)
    simulation = LandSimulation(
        start_date,
        stop_date,
        Δt,
        land;
        outdir,
        # drivers every step, so that the shortwave follows the sun
        updateat = Δt,
        diagnostics = ClimaLand.default_diagnostics(
            land,
            start_date,
            outdir;
            output_vars = OUTPUT_VARS,
            reduction_period = :monthly,
        ),
        set_ic! = (Y, p, t0, model) -> set_initial_state!(Y, ic),
    )
    @info "Segment" SCENARIO n0 n1 f0 start_date stop_date ic
    log_parameters(land, toml_dict)
    budget = CarbonBudget(simulation)
    n = n0
    last_done = n0 - 1
    year_start = time()
    try
        while n <= n1
            next_year = min(DateTime(forcing_year(n) + 1), stop_date)
            while date(simulation._integrator.t) < next_year
                ClimaLand.Simulations.step!(simulation)
                accumulate_budget!(budget, simulation, Δt)
            end
            row = budget_row!(budget, simulation)
            co2_ppm = 1e6 * co2_at(co2, DateTime(forcing_year(n), 7, 2))
            append_row(
                BUDGET_CSV,
                row;
                prefix = (;
                    nominal_year = n,
                    forcing_year = forcing_year(n),
                    co2_ppm,
                ),
            )
            save_state(simulation, n)
            row.nonfinite_nodes > 0 && error(
                "$(row.nonfinite_nodes) land nodes with non-finite carbon stocks after year $n",
            )
            push!(seconds_per_year, time() - year_start)
            @info "Year done" n row.NBP row.closure row.C_total seconds_per_year[end]
            last_done = n
            n += 1
            year_start = time()
            remaining = WALLTIME_HOURS * 3600 - (time() - JOB_START)
            if n <= n1 && (
                remaining < 1.3 * maximum(seconds_per_year) + 900 ||
                length(seconds_per_year) >= MAX_YEARS
            )
                @info "Stopping before the walltime" remaining
                break
            end
        end
    finally
        ClimaLand.Diagnostics.close_output_writers(simulation.diagnostics)
        ClimaUtilities.FileReaders.close_all_ncfiles()
    end
    return last_done
end

function main()
    context = ClimaComms.context()
    ClimaComms.init(context)
    mkpath(STATE_DIR)
    toml_dict = LP.create_toml_dict(FT)
    domain = ClimaLand.Domains.global_box_domain(
        FT;
        context,
        mask_threshold = FT(0.99),
    )
    seconds_per_year = Float64[]
    last = last_saved_year()
    n = isnothing(last) ? START_YEAR : last + 1
    n > END_YEAR && return @info "Run already complete" END_YEAR
    ic = isnothing(last) ? IC : state_path(last)
    isfile(ic) && isfile(CO2_FILE) || error(
        "set ARTIFACT_DIR to the folder of the climaland_1850_spunup_ic artifact (or IC and CO2_FILE)",
    )
    while n <= END_YEAR
        n1 = segment_end(n, END_YEAR)
        done = run_segment(n, n1, domain, toml_dict, ic, seconds_per_year)
        (done < n1 || length(seconds_per_year) >= MAX_YEARS) && break
        n = n1 + 1
        ic = state_path(n1)
        GC.gc()
    end
    @info "Job finished" last_saved_year() END_YEAR
end

main()
