# # TRENDY v15 (Global Carbon Budget 2026) simulations with ClimaLand
#
# Runs one TRENDY simulation and writes it in the TRENDY format with
# `trendy_output.jl` (in this folder):
# `<outdir>/<MODEL>/<SIM>/<MODEL>_<SIM>_<var>.nc`, the zonal NBP file and a README.
#
# Configurations
# - `A` (model name `ClimaLand-OBS`): vegetation from observations. Leaf area index
#   (LAI) is prescribed from MODIS (Yuan et al. 2011, Collection 6.1, 2000-2025);
#   the live vegetation carbon is carried by prognostic pools (sugar, leaf, wood,
#   fine roots) that follow the GPP of the prescribed LAI (`PrognosticCarbonModel`
#   around `PrescribedBiomassModel`). No biomass observations constrain the pools
#   yet. Meant for 2003-2025.
# - `B` (model name `ClimaLand`): prognostic vegetation from 1700. Not implemented.
#
# Simulations: `S2` (time-varying climate and CO2, fixed land use) is implemented.
# `S3` (with land use) is accepted but stops with an error until the removal of
# vegetation carbon by forest loss is implemented.
#
# Forcing
# - Climate: CRU-JRA v4.10 (`prescribed_forcing_crujra(...; version = "4.10")`) with
#   the default time handling (shortwave by cosine of the zenith angle, longwave
#   and precipitation as interval means, harmonic interpolation of the states)
#   and 1-hour driver updates in both the reader and the simulation.
# - CO2: the annual global means of the protocol (`prescribed_co2_annual`).
#
# Spin-up, in a separate simulation without output that starts at `--start`:
# 1. the vegetation carbon pools are set to their equilibrium under the mean
#    drivers of the first `--veg-spinup-months` months (default 12),
#    `spin_up_carbon_pools!`;
# 2. over the next `--soil-spinup-months` months (default 12), the soil organic
#    carbon is held at the SoilGrids stock of the initial conditions and the
#    decomposition rate of each soil layer is set so that its respiration
#    equals its litter input, `balance_soil_carbon!`. As in CARDAMOM, the stock
#    is observed and the turnover follows from it: the SOC has no spin-up trend.
# The TRENDY simulation then starts at `--start` from the state reached at the
# end of the spin-up (soil water, temperature, snow and carbon pools) with the
# rates of step 2.
#
# Usage, from the repository root (GPU strongly recommended):
# ```
# julia --project=.buildkite experiments/trendy/run_trendy.jl \
#     --configuration A --simulation S2 --start 2003-01-01 --stop 2026-01-01 \
#     --outdir /path/to/output
# ```
# Options (`--help` lists them all):
# - `--resolution` of the TRENDY files in degrees (1.0, the default, or 0.5);
# - `--crujra-path`: folder with `crujra_forcing_data_YYYY_0.5x0.5.nc` of CRU-JRA
#   v4.10 (default: environment variable `CRUJRA_V410_PATH`, else
#   `/net/sampo/data1/crujra/crujra_v4.10_forcing_data_artifact`);
# - `--co2-file`: annual CO2 text file (default: `TRENDY_CO2_FILE`, else
#   `/net/sampo/data1/trendy/v15/co2/global_co2_ann_1700_2025.txt`);
# - `--lai-path`: folder with `Yuan_et_al_YYYY_0.25x0.25.nc` of MODIS LAI v6.1
#   (default: `MODIS_LAI_V061_PATH`, else the `modis_lai_yuan_v061_0.25deg`
#   artifact, found by its hash in the Julia depots, honouring `Overrides.toml`;
#   this artifact is not yet in ClimaLand's `Artifacts.toml`, see
#   CliMA/ClimaArtifacts#199);
# - `--convert-only`: write the TRENDY files from existing diagnostics.
#
# Output: `<outdir>/diagnostics/` (monthly ClimaLand diagnostics) and
# `<outdir>/trendy/<MODEL>/<SIM>/` (TRENDY files).

import ClimaComms
ClimaComms.@import_required_backends
import ClimaDiagnostics
import ClimaParams as CP
using ClimaCore
using ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
import ArgParse
import Artifacts
using Dates

include(joinpath(@__DIR__, "trendy_output.jl"))

const FT = Float64

"Git tree hash of the `modis_lai_yuan_v061_0.25deg` artifact (CliMA/ClimaArtifacts#199)."
const MODIS_LAI_V061_TREE_HASH = "94a53d2e0ea64868305cd123068a71fa4ada9e11"
const DEFAULT_CRUJRA_V410_PATH = "/net/sampo/data1/crujra/crujra_v4.10_forcing_data_artifact"
const DEFAULT_CO2_FILE = "/net/sampo/data1/trendy/v15/co2/global_co2_ann_1700_2025.txt"
const MODEL_NAMES = Dict("A" => "ClimaLand-OBS", "B" => "ClimaLand")

function parse_commandline(args = ARGS)
    s = ArgParse.ArgParseSettings(
        description = "Run a TRENDY v15 simulation with ClimaLand and write the TRENDY files",
    )
    ArgParse.@add_arg_table! s begin
        "--configuration"
        help = "A (ClimaLand-OBS, prescribed LAI) or B (ClimaLand, prognostic; not implemented)"
        default = "A"
        "--simulation"
        help = "TRENDY simulation: S2 (S3 is not implemented yet)"
        default = "S2"
        "--start"
        help = "first date simulated (and written)"
        default = "2003-01-01"
        "--stop"
        help = "date at which the simulation stops (exclusive); whole years for TRENDY files"
        default = "2026-01-01"
        "--resolution"
        help = "resolution of the TRENDY files in degrees (1.0 or 0.5)"
        arg_type = Float64
        default = 1.0
        "--outdir"
        help = "output folder (diagnostics/ and trendy/ are created in it)"
        required = true
        "--crujra-path"
        help = "folder of the CRU-JRA v4.10 yearly files"
        default = get(ENV, "CRUJRA_V410_PATH", DEFAULT_CRUJRA_V410_PATH)
        "--co2-file"
        help = "annual CO2 file of the protocol"
        default = get(ENV, "TRENDY_CO2_FILE", DEFAULT_CO2_FILE)
        "--lai-path"
        help = "folder of the MODIS LAI v6.1 yearly files (default: artifact by hash)"
        default = get(ENV, "MODIS_LAI_V061_PATH", "")
        "--veg-spinup-months"
        help = "months of the vegetation-pool spin-up before the start (0: none)"
        arg_type = Int
        default = 12
        "--soil-spinup-months"
        help = "months of the soil-carbon balance after the vegetation spin-up (0: none, SOC drifts)"
        arg_type = Int
        default = 12
        "--dt"
        help = "time step in seconds"
        arg_type = Float64
        default = 900.0
        "--updateat-hours"
        help = "driver update interval in hours (the reader and the simulation use the same)"
        arg_type = Int
        default = 1
        "--convert-only"
        help = "only write the TRENDY files from <outdir>/diagnostics"
        action = :store_true
        "--allow-partial-years"
        help = "also write years that are not complete (short tests)"
        action = :store_true
    end
    return ArgParse.parse_args(args, s)
end

"""
    modis_lai_v061_dir(path)

Folder of the MODIS LAI v6.1 files: `path` if not empty, else the
`modis_lai_yuan_v061_0.25deg` artifact found by its hash (overrides honoured).
"""
function modis_lai_v061_dir(path)
    isempty(path) || return path
    dir = Artifacts.artifact_path(Base.SHA1(MODIS_LAI_V061_TREE_HASH))
    isdir(dir) || error(
        "MODIS LAI v6.1 not found at $dir. Install the modis_lai_yuan_v061_0.25deg " *
        "artifact (CliMA/ClimaArtifacts#199) or pass --lai-path / MODIS_LAI_V061_PATH.",
    )
    return dir
end

"""
    modis_lai_v061_paths(dir, start_date, stop_date)

Yearly MODIS LAI files covering `start_date` to `stop_date`, chosen as in
`ClimaLand.Artifacts.modis_lai_multiyear_paths`.
"""
function modis_lai_v061_paths(dir, start_date, stop_date)
    year0 = Dates.year(start_date)
    yearf = Dates.year(stop_date) + (Dates.month(stop_date) == 12 ? 1 : 0)
    paths = [joinpath(dir, "Yuan_et_al_$(y)_0.25x0.25.nc") for y in year0:yearf]
    for path in paths
        isfile(path) || error("MODIS LAI file $path does not exist")
    end
    return paths
end

"""
    configuration_a_model(domain, start_date, stop_date, Δt, toml_dict, options)

`LandModel` of configuration A: CRU-JRA v4.10, TRENDY CO2, prescribed MODIS LAI
v6.1 with prognostic carbon pools.
"""
function configuration_a_model(
    domain,
    start_date,
    stop_date,
    Δt,
    toml_dict,
    options,
)
    surface_space = domain.space.surface
    context = ClimaComms.context(surface_space)
    c_co2 = ClimaLand.prescribed_co2_annual(options["co2-file"], start_date)
    atmos, radiation = ClimaLand.prescribed_forcing_crujra(
        start_date,
        stop_date,
        surface_space,
        toml_dict,
        FT;
        c_co2,
        version = "4.10",
        crujra_path = options["crujra-path"],
        updateat = Dates.Hour(options["updateat-hours"]),
        context,
    )
    lai_paths = modis_lai_v061_paths(
        modis_lai_v061_dir(options["lai-path"]),
        start_date,
        stop_date,
    )
    LAI = ClimaLand.Canopy.prescribed_lai_modis(
        surface_space,
        start_date,
        stop_date;
        modis_lai_ncdata_path = lai_paths,
    )
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    biomass = ClimaLand.Canopy.PrognosticCarbonModel{FT}(
        ClimaLand.Canopy.PrescribedBiomassModel{FT}(
            surface_domain,
            LAI,
            ClimaLand.Canopy.modis_max_lai(surface_space),
            toml_dict,
        ),
        toml_dict;
        tree_share = ClimaLand.Canopy.clm_tree_share(surface_space),
    )
    return LandModel{FT}(
        (; atmos, radiation),
        LAI,
        toml_dict,
        domain,
        Δt;
        prognostic_land_components = (:canopy, :lake, :snow, :soil, :soilco2),
        biomass,
    ),
    lai_paths
end

"ClimaLand diagnostics needed by the TRENDY writer, plus a few for checks."
function trendy_diagnostic_names(model)
    possible = ClimaLand.Diagnostics.get_possible_diagnostics(model)
    wanted = unique(
        vcat(
            [v.inputs for v in trendy_variables()]...,
            ["er", "swe", "soc_int"],
        ),
    )
    missing_names = setdiff(wanted, possible)
    isempty(missing_names) ||
        @warn "Diagnostics not available in this model" missing_names
    return filter(in(possible), wanted)
end

function git_commit()
    try
        return readchomp(`git -C $(pkgdir(ClimaLand)) rev-parse HEAD`)
    catch
        return "unknown"
    end
end

function main(options)
    configuration = options["configuration"]
    simulation_name = options["simulation"]
    configuration == "B" && error(
        "Configuration B (ClimaLand, prognostic vegetation from 1700) is not implemented yet",
    )
    configuration == "A" ||
        error("Unknown configuration $configuration (A or B)")
    simulation_name == "S3" && error(
        "S3 is not implemented yet for configuration A: forest-loss removal not implemented yet",
    )
    simulation_name == "S2" || error(
        "Simulation $simulation_name is not implemented for configuration A (S2 only; S3 to come)",
    )
    model_name = MODEL_NAMES[configuration]
    start_date = DateTime(options["start"])
    stop_date = DateTime(options["stop"])
    outdir = options["outdir"]
    diagnostics_dir = joinpath(outdir, "diagnostics")
    trendy_dir = joinpath(outdir, "trendy")
    resolution = options["resolution"]
    domain_depth = 15.0

    if !options["convert-only"]
        context = ClimaComms.context()
        ClimaComms.init(context)
        Δt = options["dt"]
        updateat = Dates.Hour(options["updateat-hours"])
        domain = ClimaLand.Domains.global_box_domain(
            FT;
            context,
            mask_threshold = FT(0.99),
            depth = domain_depth,
        )
        toml_dict = LP.create_toml_dict(FT)
        spinup_months = options["veg-spinup-months"]
        soil_months = options["soil-spinup-months"]
        soil_months > 0 &&
            spinup_months == 0 &&
            error("The soil-carbon balance needs the vegetation spin-up first")
        veg_end = start_date + Month(spinup_months)
        spinup_end = veg_end + Month(soil_months)
        # The drivers cover both the spin-up and the simulation
        model, lai_paths = configuration_a_model(
            domain,
            start_date,
            max(stop_date, spinup_end),
            Δt,
            toml_dict,
            options,
        )
        @info "TRENDY run" model_name simulation_name start_date stop_date resolution spinup_months updateat options["crujra-path"] options["co2-file"] lai_paths soil_months
        state = nothing
        ln_rate_factor = nothing
        if spinup_months > 0
            spinup = LandSimulation(
                start_date,
                spinup_end,
                Δt,
                model;
                outdir,
                diagnostics = (),
                updateat,
            )
            # SoilGrids SOC of the initial conditions, the stock held by step 2
            SOC_initial = copy(spinup._integrator.u.soilco2.SOC)
            t0 = time()
            ClimaLand.Simulations.spin_up_carbon_pools!(spinup, veg_end)
            @info "Vegetation-pool spin-up done" veg_end wall_minutes =
                round((time() - t0) / 60; digits = 1)
            if soil_months > 0
                t0 = time()
                spinup._integrator.u.soilco2.SOC .= SOC_initial
                ClimaLand.Simulations.balance_soil_carbon!(spinup, spinup_end)
                ln_rate_factor =
                    copy(spinup._integrator.p.soilco2.ln_rate_factor)
                # Distribution of the factor over the soil layers with SOC
                # (ocean points of the mask and unchanged layers have 0)
                lf = vec(Array(parent(ln_rate_factor)))
                soc = vec(Array(parent(SOC_initial)))
                set = @. isfinite(soc) & (soc > 0) & (lf != 0)
                f = sort(exp.(lf[set]))
                q(x) = round(
                    f[clamp(round(Int, x * length(f)), 1, length(f))];
                    sigdigits = 3,
                )
                @info "Soil-carbon balance done" spinup_end wall_minutes =
                    round((time() - t0) / 60; digits = 1) n_layers = length(f) factor_quantiles_1_10_50_90_99 =
                    q.((0.01, 0.1, 0.5, 0.9, 0.99)) share_below_0_01 =
                    round(count(<(0.01), f) / length(f); digits = 3) share_above_100 =
                    round(count(>(100), f) / length(f); digits = 3)
            end
            state = copy(spinup._integrator.u)
            spinup = nothing
            GC.gc()
        end
        set_ic! =
            isnothing(state) ?
            ClimaLand.Simulations.make_set_initial_state_from_file(
                ClimaLand.Artifacts.rosetta_spunup_ic_path(; context),
                model,
            ) :
            function (Y, p, t0, model)
                Y .= state
                isnothing(ln_rate_factor) ||
                    (p.soilco2.ln_rate_factor .= ln_rate_factor)
                return nothing
            end
        mkpath(diagnostics_dir)
        space = domain.space.subsurface
        output_writer = ClimaDiagnostics.Writers.NetCDFWriter(
            space,
            diagnostics_dir;
            start_date,
            horizontal_pts = trendy_grid(resolution),
        )
        diagnostics = ClimaLand.default_diagnostics(
            model,
            start_date,
            diagnostics_dir;
            output_writer,
            output_vars = trendy_diagnostic_names(model),
            reduction_period = :monthly,
        )
        simulation = LandSimulation(
            start_date,
            stop_date,
            Δt,
            model;
            outdir = diagnostics_dir,
            set_ic!,
            diagnostics,
            updateat,
        )
        t0 = time()
        solve!(simulation)
        wall = time() - t0
        years = Dates.value(stop_date - start_date) / (1000 * 86400 * 365.25)
        @info "Simulation done" wall_hours = round(wall / 3600; digits = 2) wall_minutes_per_simulated_year =
            round(wall / 60 / years; digits = 1)
    end

    written = write_trendy_output(
        diagnostics_dir,
        trendy_dir;
        model = model_name,
        simulation = simulation_name,
        resolution,
        commit = git_commit(),
        soil_depth = domain_depth,
        allow_partial_years = options["allow-partial-years"],
    )
    @info "TRENDY files written" trendy_dir length(written)
    return written
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(parse_commandline())
end
