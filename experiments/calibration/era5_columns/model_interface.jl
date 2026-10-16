# Forward model and observation map of the ERA5 column-ensemble calibration:
# the land model of the default long run (P-model, prescribed MODIS LAI) on an
# ensemble of columns, forced by ERA5.

import ClimaComms
ClimaComms.@import_required_backends
import ClimaCalibrate
import ClimaCore
import ClimaDiagnostics
import ClimaParams as CP
import ClimaUtilities.TimeManager: date
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
import Dates
import EnsembleKalmanProcesses as EKP
import JLD2
import TOML

include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "observations.jl"))
include(joinpath(@__DIR__, "columns.jl"))

"""Name of the file holding the seasonal means simulated by a run."""
const SEASONAL_MEANS_FILE = "seasonal_means.jld2"

"""Name of the file, in the calibration directory, holding the observations
and which of their entries are calibrated."""
const OBSERVATIONS_FILE = "observations.jld2"

"""Converts a simulated diagnostic to the units of its observations: GPP from
mol CO2 m^-2 s^-1 to g C m^-2 day^-1. The other diagnostics already match."""
const SIM_TO_OBS_UNITS = Dict("gpp" => x -> x * 86400 * 12.011)

"""
    ERA5ColumnsInterface <: ClimaCalibrate.AbstractModelInterface

Model interface of the ERA5 column-ensemble calibration.
"""
struct ERA5ColumnsInterface <: ClimaCalibrate.AbstractModelInterface
    "Directory of the iterations, the observations, and the results"
    output_dir::String
    "Whether the forcing is the 8° ERA5 of 2008, repeated every year, instead
    of the hourly 1° ERA5 of the simulated years"
    lowres_forcing::Bool
end

ClimaCalibrate.model_interface_filepath(::ERA5ColumnsInterface) = @__FILE__

"""
    columns_path(interface::ERA5ColumnsInterface)

Return the path to the CSV file of the simulated columns.
"""
columns_path(interface::ERA5ColumnsInterface) =
    joinpath(interface.output_dir, "results", "columns.csv")

"""
    setup_simulation(
        toml_dict,
        longlat;
        lowres_forcing = false,
        start_date = SPINUP_START,
        stop_date = STOP_DATE,
    )

Return the simulation of the default land model on columns at `longlat`, forced
by ERA5, and the `DictWriter` collecting the monthly means of `SHORT_NAMES`.
"""
function setup_simulation(
    toml_dict,
    longlat;
    lowres_forcing = false,
    start_date = SPINUP_START,
    stop_date = STOP_DATE,
)
    context = ClimaComms.context()
    domain = ClimaLand.Domains.ColumnEnsemble(;
        zlim = (FT(-15), FT(0)),
        nelements = 15,
        dz_tuple = (FT(3), FT(0.05)),
        longlat = [FT.(ll) for ll in longlat],
    )
    surface_space = domain.space.surface
    check_column_order(surface_space, longlat)

    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        start_date,
        stop_date,
        surface_space,
        toml_dict,
        FT;
        max_wind_speed = 25.0,
        use_lowres_forcing = lowres_forcing,
        context,
    )
    LAI = ClimaLand.Canopy.prescribed_lai_modis(
        surface_space,
        start_date,
        stop_date,
    )
    model = ClimaLand.LandModel{FT}(
        (; atmos, radiation),
        LAI,
        toml_dict,
        domain,
        Δt;
        prognostic_land_components = (:canopy, :lake, :snow, :soil, :soilco2),
    )

    writer = ClimaDiagnostics.Writers.DictWriter()
    diagnostics = ClimaLand.default_diagnostics(
        model,
        start_date,
        nothing;
        output_writer = writer,
        output_vars = SHORT_NAMES,
        reduction_period = :monthly,
        reduction_type = :average,
    )
    simulation = LandSimulation(start_date, stop_date, Δt, model; diagnostics)
    return simulation, writer
end

"""
    check_column_order(surface_space, longlat)

Error unless the points of `surface_space` are at `longlat`, in order, since
the outputs are matched to the observations by position.
"""
function check_column_order(surface_space, longlat)
    coordinates = ClimaCore.Fields.coordinate_field(surface_space)
    longs = vec(Array(parent(coordinates.long)))
    lats = vec(Array(parent(coordinates.lat)))
    (
        longs ≈ [FT(lon) for (lon, _) in longlat] && lats ≈ [FT(lat) for (_, lat) in longlat]
    ) || error("The columns are not in the order of `longlat`")
    return nothing
end

"""
    monthly_means(writer, short_name)

Return the months (first days) and the `months × columns` matrix of the
monthly means of `short_name` saved in `writer`, in observation units.
"""
function monthly_means(writer, short_name)
    saved = writer["$(short_name)_1M_average"]
    times = sort!(collect(keys(saved)))
    # A monthly mean is saved at the start of the following month
    months = [Dates.Date(date(t)) - Dates.Month(1) for t in times]
    values =
        reduce(vcat, [permutedims(vec(Array(parent(saved[t])))) for t in times])
    return months, get(SIM_TO_OBS_UNITS, short_name, identity).(values)
end

"""
    simulated_seasonal_means(writer)

Return a dictionary mapping each of `SHORT_NAMES` to the `seasons × columns`
matrix of its seasonal means over the calibration period.
"""
function simulated_seasonal_means(writer)
    seasons = season_starts(CALIBRATION_START, STOP_DATE)
    return Dict(map(SHORT_NAMES) do short_name
        months, monthly = monthly_means(writer, short_name)
        short_name => seasonal_means(monthly, months, seasons)
    end)
end

"""
    write_parameter_toml(path, parameters; header = "")

Write the `parameters` (a dictionary from names to values) to the TOML file at
`path`, in the format of the ClimaLand parameter files, preceded by the comment
lines in `header`.
"""
function write_parameter_toml(path, parameters; header = "")
    open(path, "w") do io
        print(io, header)
        for name in sort!(collect(keys(parameters)))
            println(io, "\n[\"$name\"]")
            println(io, "value = $(parameters[name])")
            println(io, "type = \"float\"")
        end
    end
    return nothing
end

"""
    model_parameter_file(parameter_file)

Return the path to a TOML file holding the model parameters that the
calibrated parameters in `parameter_file` set (see `model_parameters`): a file
written next to it if it holds `STEM_REFLECTANCE_SCALE`, or `parameter_file`
itself.
"""
function model_parameter_file(parameter_file)
    entries = TOML.parsefile(parameter_file)
    haskey(entries, STEM_REFLECTANCE_SCALE) || return parameter_file
    calibrated = Dict(name => entry["value"] for (name, entry) in entries)
    path = joinpath(dirname(parameter_file), "model_parameters.toml")
    write_parameter_toml(
        path,
        model_parameters(calibrated, LP.create_toml_dict(FT)),
    )
    return path
end

"""
    run_columns(interface::ERA5ColumnsInterface, parameter_files, path)

Run the columns of `interface` with the default parameters overridden by the
calibrated parameters in the TOML files `parameter_files`, and save the
simulated seasonal means to `path`.
"""
function run_columns(interface::ERA5ColumnsInterface, parameter_files, path)
    override_files = model_parameter_file.(parameter_files)
    toml_dict = LP.create_toml_dict(FT; override_files)
    simulation, writer = setup_simulation(
        toml_dict,
        read_columns(columns_path(interface)).longlat;
        interface.lowres_forcing,
    )
    for file in override_files
        CP.check_override_parameter_usage(
            toml_dict,
            keys(TOML.parsefile(file)),
            true,
        )
    end
    solve!(simulation)
    mkpath(dirname(path))
    JLD2.save_object(path, simulated_seasonal_means(writer))
    return nothing
end

function ClimaCalibrate.forward_model(
    interface::ERA5ColumnsInterface,
    iteration,
    member,
)
    (; output_dir) = interface
    member_path =
        ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, member)
    run_columns(
        interface,
        [ClimaCalibrate.parameter_path(output_dir, iteration, member)],
        joinpath(member_path, SEASONAL_MEANS_FILE),
    )
    return nothing
end

"""
    flatten(seasonal, masks)

Stack the calibrated entries (`masks`) of the seasonal means of each of
`SHORT_NAMES` into a vector ordered as the observations.
"""
flatten(seasonal, masks) =
    reduce(vcat, [vec(seasonal[name])[masks[name]] for name in SHORT_NAMES])

function ClimaCalibrate.observation_map(
    interface::ERA5ColumnsInterface,
    iteration,
)
    (; output_dir) = interface
    masks = JLD2.load(joinpath(output_dir, OBSERVATIONS_FILE), "masks")
    ekp = ClimaCalibrate.load_ekp_struct(output_dir, iteration)
    n_observations = sum(sum, values(masks))
    g_columns = map(1:EKP.get_N_ens(ekp)) do member
        path = joinpath(
            ClimaCalibrate.path_to_ensemble_member(
                output_dir,
                iteration,
                member,
            ),
            SEASONAL_MEANS_FILE,
        )
        if isfile(path)
            flatten(JLD2.load_object(path), masks)
        else
            @error "Member $member has no output, filling it with NaNs"
            fill(NaN, n_observations)
        end
    end
    return reduce(hcat, g_columns)
end
