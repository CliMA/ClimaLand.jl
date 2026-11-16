# TRENDY-style global runs of the LandModel with prognostic LAI and carbon pools.
#
# A run covers the nominal years NOMINAL_START to NOMINAL_END (inclusive). Each nominal
# year n is forced by the forcing year f(n): with RECYCLE=1, or for n < 1901, the
# 1901-1920 cycle, f(n) = 1901 + mod(n - CYCLE_ORIGIN, 20); otherwise f(n) = n.
# Consecutive nominal years with consecutive forcing years form a segment, run as one
# LandSimulation from the state at the end of the previous year. CO2 follows the nominal
# year (CO2_CSV), or stays at CO2_FIXED_YEAR's annual mean.
#
# At the end of each year the state is saved to states/state_<n>.hdf5 and a row of
# global carbon fluxes and stocks is appended to annual_budget.csv. A rerun resumes after
# the last saved year, and a job stops at a year boundary before WALLTIME_HOURS.
#
# Options (environment variables), with defaults:
#   NOMINAL_START, NOMINAL_END     required
#   FORCING=crujra                 or era5
#   RECYCLE=0, CYCLE_ORIGIN=1850
#   CO2_CSV=co2_annual_1700_2025.csv (this folder), CO2_FIXED_YEAR (unset)
#   IC_FILE (unset: the default initial conditions)
#   OVERRIDE_TOML (unset), WALLTIME_HOURS=12, MAX_YEARS (unset: no limit per job)
#   OUTPUT_VARS (comma-separated short names; default below)
#   SPINUP (unset, `pools` or `pools_soc`): accumulate the spin-up sums over the whole run
#     (spinup_tools.jl; saved with each state as states/accum_<n>.hdf5) and, after
#     NOMINAL_END, write states/spun_<NOMINAL_END>.hdf5 with the vegetation pools (and
#     with `pools_soc` the SOC) at equilibrium, and spinup_report.txt; SOC_CAP=100
#   SOC_ALTS (unset): alternative DAMM curve sets for the SOC sums, e.g.
#     "ea=60000;kmo2=0.3;ea=45000,kmo2=0.1;doa=3" (spinup_tools.jl)
#   APPLY_SOC=1: with SPINUP=pools_soc, set the SOC to the solve at the end (0: keep it,
#     for an offline choice of the DAMM parameters from the sums)
#   REPAIR_IC=0: with 1, columns of the IC file holding a non-finite value take the
#     values of the default initial conditions
#   The accumulators of all but the last two years are deleted.
#   NAN_MAX=0: stop with an error after a year that ends with more land nodes than this
#     where a carbon stock is not finite
#   PROBE_NODE, PROBE_FROM (unset): record every step from the date PROBE_FROM all real
#     fields of Y and p at the surface node PROBE_NODE (linear index in parent arrays;
#     probe.jl), and write probe_<n0>.csv at the end of the segment, or with an error at
#     the first non-finite value of Y there

import ClimaComms
ClimaComms.@import_required_backends
import ClimaCore
import ClimaUtilities
import ClimaUtilities.TimeManager: ITime, date
import ClimaParams as CP
using ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation
using Dates

include(joinpath(@__DIR__, "co2_driver.jl"))
include(joinpath(@__DIR__, "carbon_budget.jl"))
include(joinpath(@__DIR__, "spinup_tools.jl"))
include(joinpath(@__DIR__, "state_io.jl"))
include(joinpath(@__DIR__, "probe.jl"))

const FT = Float64
const Δt = 900.0
const DEFAULT_OUTPUT_VARS = [
    "gpp", "nee", "er", "hr", "ra", "lai", "cveg", "soc_int", "ftr",
    "lhf", "shf", "swu", "lwu", "et", "tair", "precip", "swd",
]

env(name, default) = get(ENV, name, default)
env_int(name, default) = parse(Int, env(name, string(default)))

const NOMINAL_START = env_int("NOMINAL_START", 0)
const NOMINAL_END = env_int("NOMINAL_END", 0)
NOMINAL_START > 0 && NOMINAL_END >= NOMINAL_START ||
    error("set NOMINAL_START and NOMINAL_END")
const FORCING = env("FORCING", "crujra")
const RECYCLE = env("RECYCLE", "0") == "1"
const CYCLE_ORIGIN = env_int("CYCLE_ORIGIN", 1850)
const CO2_CSV = env("CO2_CSV", joinpath(@__DIR__, "co2_annual_1700_2025.csv"))
const CO2_FIXED_YEAR = haskey(ENV, "CO2_FIXED_YEAR") ? parse(Int, ENV["CO2_FIXED_YEAR"]) : nothing
const IC_FILE = env("IC_FILE", "")
const OVERRIDE_TOML = env("OVERRIDE_TOML", "")
const WALLTIME_HOURS = parse(Float64, env("WALLTIME_HOURS", "12"))
const MAX_YEARS = haskey(ENV, "MAX_YEARS") ? parse(Int, ENV["MAX_YEARS"]) : typemax(Int)
const SPINUP = env("SPINUP", "")
SPINUP in ("", "pools", "pools_soc") || error("unknown SPINUP $SPINUP")
const SOC_CAP = parse(Float64, env("SOC_CAP", "100"))
const NAN_MAX = env_int("NAN_MAX", 0)
const SOC_ALTS = parse_alts(env("SOC_ALTS", ""))
const APPLY_SOC = env("APPLY_SOC", "1") == "1"
const REPAIR_IC = env("REPAIR_IC", "0") == "1"
const PROBE_NODE = haskey(ENV, "PROBE_NODE") ? parse(Int, ENV["PROBE_NODE"]) : nothing
const PROBE_FROM = DateTime(env("PROBE_FROM", "1700-01-01"))
const OUTPUT_VARS = haskey(ENV, "OUTPUT_VARS") ? split(ENV["OUTPUT_VARS"], ',') .|> String : DEFAULT_OUTPUT_VARS
const RUN_DIR = pwd()
const STATE_DIR = joinpath(RUN_DIR, "states")
const BUDGET_CSV = joinpath(RUN_DIR, "annual_budget.csv")
const JOB_START = time()

forcing_year(n) = (RECYCLE || n < 1901) ? 1901 + mod(n - CYCLE_ORIGIN, 20) : n

"""Consecutive nominal years from `n0` to `n1` with consecutive forcing years."""
function segment_end(n0, n1)
    n = n0
    while n < n1 && forcing_year(n + 1) == forcing_year(n) + 1
        n += 1
    end
    return n
end

state_path(n::Integer) = joinpath(STATE_DIR, "state_$(n).hdf5")
state_path(name::AbstractString) = joinpath(STATE_DIR, "$(name).hdf5")
accum_path(n) = joinpath(STATE_DIR, "accum_$(n).hdf5")

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
    ClimaLand.save_checkpoint(simulation._integrator.u, simulation._integrator.t, tmp; model = simulation.model)
    mv(only(readdir(tmp; join = true)), state_path(n); force = true)
    rm(tmp; recursive = true)
end

function make_land(start_date, stop_date, domain, toml_dict, co2)
    surface_space = domain.space.surface
    context = ClimaComms.context(surface_space)
    c_co2 = co2_time_varying_input(co2)
    forcing = if FORCING == "crujra"
        ClimaLand.prescribed_forcing_crujra(start_date, stop_date, surface_space, toml_dict, FT; c_co2, context, fill_missing = true, zenith_scaled_shortwave = true)
    elseif FORCING == "era5"
        ClimaLand.prescribed_forcing_era5(start_date, stop_date, surface_space, toml_dict, FT; max_wind_speed = 25.0, c_co2, context)
    else
        error("unknown FORCING $FORCING")
    end
    biomass = ClimaLand.Canopy.PrognosticCarbonModel{FT}(
        ClimaLand.Canopy.ZhouOptimalLAIModel{FT}(ClimaLand.Domains.obtain_surface_domain(domain), toml_dict),
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

"""
    repaired_ic(ic_file, land)

A `set_ic!` that reads the state of `ic_file`, then gives every surface node where a
field of it is not finite the values of the default initial conditions in all fields.
"""
function repaired_ic(ic_file, land)
    default_ic! = ClimaLand.Simulations.make_set_initial_state_from_file(
        ClimaLand.Artifacts.rosetta_spunup_ic_path(; context = ClimaComms.context(land)),
        land,
    )
    function set_ic!(Y, p, t0, model)
        default_ic!(Y, p, t0, model)
        Y0 = copy(Y)
        read_into!(Y, ic_file, "Y")
        leaves(x) = x isa ClimaCore.Fields.Field ? [x] : reduce(vcat, (leaves(getproperty(x, k)) for k in propertynames(x)); init = Any[])
        columns(f) = (a = Array(parent(f)); reshape(a, :, size(a, ndims(a))))
        bad = falses(0)
        for f in leaves(Y)
            c = columns(f)
            isempty(bad) && (bad = falses(size(c, 2)))
            bad .|= vec(.!all(isfinite, c; dims = 1))
        end
        mask = ClimaCore.Spaces.get_mask(land.soil.domain.space.surface)
        mask isa ClimaCore.DataLayouts.NoMask || (bad .&= vec(Array(parent(mask.is_active))))
        for (f, f0) in zip(leaves(Y), leaves(Y0))
            c, c0 = columns(f), columns(f0)
            c[:, bad] .= c0[:, bad]
            copyto!(parent(f), reshape(c, size(parent(f))))
        end
        @info "IC repaired from the default initial conditions" count(bad) nodes = findall(bad)
    end
    return set_ic!
end

"""
    log_parameters(land, toml_dict)

Writes parameters_used.txt, the soil CO2 (DAMM) and carbon pool parameters as built,
and parameters.toml, the TOML parameters the model read, once per run directory.
"""
function log_parameters(land, toml_dict)
    path = joinpath(RUN_DIR, "parameters_used.txt")
    isfile(path) && return nothing
    lines = String["override: " * (isempty(OVERRIDE_TOML) ? "none" : OVERRIDE_TOML)]
    for (label, params) in (("soilco2", land.soilco2.parameters), ("carbon pools", land.canopy.biomass.parameters))
        for name in fieldnames(typeof(params))
            value = getfield(params, name)
            value isa Number && push!(lines, "$label.$name = $value")
        end
    end
    write(path, join(lines, "\n") * "\n")
    CP.log_parameter_information(toml_dict, joinpath(RUN_DIR, "parameters.toml"))
    @info "Parameters" lines
end

"""
    finish_spinup(spin, simulation, budget, n)

Sets the vegetation pools (and with `SPINUP=pools_soc` the SOC) to their equilibrium
under the sums of the cycle, saves the state as states/spun_<n>.hdf5, and writes
spinup_report.txt.
"""
function finish_spinup(spin, simulation, budget, n)
    before = global_stocks(budget, simulation)
    apply_pool_equilibrium!(spin, simulation)
    lines = ["spin-up applied after nominal year $n ($(SPINUP)); cycle of $(spin.pool_seconds / 86400 / 365) years"]
    if SPINUP == "pools_soc" && !APPLY_SOC
        push!(lines, "SOC kept (APPLY_SOC=0); the sums are in accum_$(n).hdf5")
    elseif SPINUP == "pools_soc"
        soc, capped = solve_soc(spin, simulation; cap = SOC_CAP)
        Y = simulation._integrator.u
        capped_soc = @. soc * capped
        column_integral_definite!(budget.scratch, capped_soc)
        capped_Pg = land_integral(budget, budget.scratch) * 1e-12
        column_integral_definite!(budget.scratch, capped)
        capped_column = @. ifelse(budget.scratch > 0, 1.0, 0.0)
        capped_area = land_integral(budget, capped_column) / budget.land_area
        Y.soilco2.SOC .= soc
        push!(lines, "SOC cap $(SOC_CAP) kg C m-3: carbon in capped layers $(capped_Pg) Pg C; land share of columns with a capped layer $(capped_area)")
    end
    after = global_stocks(budget, simulation)
    for k in keys(after)
        push!(lines, "$k: $(getproperty(before, k)) -> $(getproperty(after, k)) Pg C")
    end
    write(joinpath(RUN_DIR, "spinup_report.txt"), join(lines, "\n") * "\n")
    @info "Spin-up applied" lines
    budget.last_stocks = after
end

"""
    finish_only(domain, toml_dict)

Applies the spin-up to the state saved after NOMINAL_END, with the sums saved with it,
for a run whose years are done but whose spin-up was not applied.
"""
function finish_only(domain, toml_dict)
    n = NOMINAL_END
    f = forcing_year(n)
    start_date = DateTime(f)
    co2 = AnnualCO2Series(CO2_CSV; fixed_year = something(CO2_FIXED_YEAR, 1850))
    land = make_land(start_date, DateTime(f + 1), domain, toml_dict, co2)
    simulation = LandSimulation(
        start_date,
        DateTime(f + 1),
        Δt,
        land;
        outdir = mktempdir(RUN_DIR),
        diagnostics = (),
        set_ic! = (Y, p, t0, model) -> read_into!(Y, state_path(n), "Y"),
    )
    context = ClimaComms.context(simulation.model)
    spin = SpinupAccumulator(simulation; soc = SPINUP == "pools_soc", alts = SOC_ALTS)
    load_spinup!(spin, accum_path(n), context)
    budget = CarbonBudget(simulation)
    finish_spinup(spin, simulation, budget, n)
    save_state(simulation, "spun_$(n)")
end

"""Runs the nominal years n0..n1 of one segment; returns the last year completed."""
function run_segment(n0, n1, domain, toml_dict, ic_file, seconds_per_year, spin_ref)
    f0 = forcing_year(n0)
    start_date = DateTime(f0)
    stop_date = DateTime(f0 + n1 - n0 + 1)
    co2 = isnothing(CO2_FIXED_YEAR) ? AnnualCO2Series(CO2_CSV; year_offset = n0 - f0) :
          AnnualCO2Series(CO2_CSV; fixed_year = CO2_FIXED_YEAR)
    land = make_land(start_date, stop_date, domain, toml_dict, co2)
    outdir = joinpath(RUN_DIR, "output", "seg_$(n0)_f$(f0)")
    mkpath(outdir)
    # drivers every step, so that the shortwave follows the sun
    kwargs = (;
        outdir,
        updateat = Δt,
        diagnostics = ClimaLand.default_diagnostics(land, start_date, outdir; output_vars = OUTPUT_VARS, reduction_period = :monthly),
    )
    if !isempty(ic_file)
        set_ic! = REPAIR_IC ? repaired_ic(ic_file, land) : (Y, p, t0, model) -> read_into!(Y, ic_file, "Y")
        kwargs = merge(kwargs, (; set_ic!))
    end
    simulation = LandSimulation(start_date, stop_date, Δt, land; kwargs...)
    @info "Segment" n0 n1 f0 start_date stop_date ic_file co2_start_ppm = 1e6 * co2(simulation._integrator.t)
    log_parameters(land, toml_dict)
    probe = isnothing(PROBE_NODE) ? nothing : NodeProbe(PROBE_NODE, PROBE_FROM)
    budget = CarbonBudget(simulation)
    context = ClimaComms.context(simulation.model)
    if !isempty(SPINUP) && isnothing(spin_ref[])
        spin_ref[] = SpinupAccumulator(simulation; soc = SPINUP == "pools_soc", alts = SOC_ALTS)
        isfile(accum_path(n0 - 1)) && load_spinup!(spin_ref[], accum_path(n0 - 1), context)
    end
    spin = spin_ref[]
    n = n0
    last_done = n0 - 1
    year_start = time()
    try
        while n <= n1
            next_year = DateTime(forcing_year(n) + 1)
            while date(simulation._integrator.t) < next_year
                ClimaLand.Simulations.step!(simulation)
                accumulate_budget!(budget, simulation, Δt)
                isnothing(spin) || accumulate_spinup!(spin, simulation, Δt)
                if !isnothing(probe) && probe_step!(probe, simulation)
                    write_probe(probe, joinpath(RUN_DIR, "probe_$(n0).csv"))
                    error("non-finite state at the probe node $(PROBE_NODE) on $(date(simulation._integrator.t))")
                end
            end
            row = budget_row!(budget, simulation)
            co2_ppm = 1e6 * co2_at(co2, DateTime(forcing_year(n), 7, 2))
            append_row(BUDGET_CSV, row; prefix = (; nominal_year = n, forcing_year = forcing_year(n), co2_ppm))
            save_state(simulation, n)
            row.nonfinite_nodes > NAN_MAX &&
                error("$(row.nonfinite_nodes) land nodes with non-finite carbon stocks after year $n (NAN_MAX = $NAN_MAX)")
            if !isnothing(spin)
                save_spinup(spin, accum_path(n), context)
                isfile(accum_path(n - 2)) && rm(accum_path(n - 2))
            end
            if !isnothing(spin) && n == NOMINAL_END
                finish_spinup(spin, simulation, budget, n)
                save_state(simulation, "spun_$(n)")
            end
            elapsed_year = time() - year_start
            push!(seconds_per_year, elapsed_year)
            @info "Year done" n row.NBP row.closure row.C_total elapsed_year
            last_done = n
            n += 1
            year_start = time()
            years_this_job = length(seconds_per_year)
            remaining = WALLTIME_HOURS * 3600 - (time() - JOB_START)
            if n <= n1 && (remaining < 1.3 * maximum(seconds_per_year) + 900 || years_this_job >= MAX_YEARS)
                @info "Stopping before the walltime" remaining years_this_job
                break
            end
        end
    finally
        isnothing(probe) || write_probe(probe, joinpath(RUN_DIR, "probe_$(n0).csv"))
        ClimaLand.Diagnostics.close_output_writers(simulation.diagnostics)
        ClimaUtilities.FileReaders.close_all_ncfiles()
    end
    return last_done
end

function main()
    context = ClimaComms.context()
    ClimaComms.init(context)
    mkpath(STATE_DIR)
    toml_dict = isempty(OVERRIDE_TOML) ? LP.create_toml_dict(FT) :
                LP.create_toml_dict(FT; override_files = [OVERRIDE_TOML])
    domain = ClimaLand.Domains.global_box_domain(FT; context, mask_threshold = FT(0.99))
    seconds_per_year = Float64[]
    spin_ref = Ref{Any}(nothing)
    last = last_saved_year()
    n = isnothing(last) ? NOMINAL_START : last + 1
    if n > NOMINAL_END
        if !isempty(SPINUP) && !isfile(state_path("spun_$(NOMINAL_END)"))
            finish_only(domain, toml_dict)
        else
            @info "Run already complete" NOMINAL_END
        end
        return
    end
    ic_file = isnothing(last) ? IC_FILE : state_path(last)
    while n <= NOMINAL_END
        n1 = segment_end(n, NOMINAL_END)
        done = run_segment(n, n1, domain, toml_dict, ic_file, seconds_per_year, spin_ref)
        done < n1 && break
        n = n1 + 1
        ic_file = state_path(n1)
        GC.gc()
    end
    @info "Job finished" last_saved_year() NOMINAL_END
end

main()
