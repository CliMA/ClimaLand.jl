# Forward model and observation map of the FLUXNET site calibration: the land
# model of the default long run (P-model, prescribed MODIS LAI) on a
# `ColumnEnsemble` with one column per site, forced by the site observations.

import ClimaComms
ClimaComms.@import_required_backends
import ClimaCalibrate
import ClimaCore
import ClimaDiagnostics
import ClimaParams as CP
import ClimaUtilities.TimeManager: date
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
import Dates
import EnsembleKalmanProcesses as EKP
import JLD2
import TOML

include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "sites.jl"))

"""Name of the file holding the monthly means simulated by a run."""
const MONTHLY_MEANS_FILE = "monthly_means.jld2"

"""Name of the file, in the calibration directory, holding the observations
and which of their entries are calibrated."""
const OBSERVATIONS_FILE = "observations.jld2"

"""Name of the file, in the calibration directory, holding the forcing and the
observed monthly means of the sites (see `build_site_data`)."""
const SITE_DATA_FILE = "site_data.jld2"

"""Converts a simulated diagnostic to the units of its observations: GPP from
mol CO2 m^-2 s^-1 to g C m^-2 day^-1. The other diagnostics already match."""
const SIM_TO_OBS_UNITS = Dict("gpp" => x -> x * 86400 * 12.011)

"""
    FluxnetSitesInterface <: ClimaCalibrate.AbstractModelInterface

Model interface of the FLUXNET site calibration.
"""
struct FluxnetSitesInterface <: ClimaCalibrate.AbstractModelInterface
    "Directory of the iterations, the site data, the observations, and the results"
    output_dir::String
end

ClimaCalibrate.model_interface_filepath(::FluxnetSitesInterface) = @__FILE__

"""
    sites_path(interface::FluxnetSitesInterface)

Return the path to the CSV file of the simulated sites.
"""
sites_path(interface::FluxnetSitesInterface) =
    joinpath(interface.output_dir, "results", "sites.csv")

"""
    site_data_path(interface::FluxnetSitesInterface)

Return the path to the JLD2 file of the forcing and observations of the sites.
"""
site_data_path(interface::FluxnetSitesInterface) =
    joinpath(interface.output_dir, SITE_DATA_FILE)

"""
    column_field(values, surface_space)

Return a field on `surface_space` with one of `values` per column, in the order
of `ClimaCore.Fields.field2array`.
"""
function column_field(values, surface_space)
    field = ClimaCore.Fields.zeros(surface_space)
    device_values =
        ClimaComms.array_type(ClimaComms.device(surface_space))(FT.(values))
    vec(ClimaCore.Fields.field2array(field)) .= device_values
    return field
end

"""
    site_forcing(times, forcing, tower_height, surface_space, toml_dict)

Return the `PrescribedAtmosphere` and `PrescribedRadiativeFluxes` of the
columns of `surface_space`, one per site, from the `forcing` matrices (one row
per column) at `times` (s since `SPINUP_START`) saved by `build_site_data`.

The forcing is at the height of the tower above the ground, raised to
`MIN_CANOPY_CLEARANCE` above the CLM canopy height of the column if needed.
"""
function site_forcing(times, forcing, tower_height, surface_space, toml_dict)
    column_input(name) =
        TimeVaryingInput(times, getproperty(forcing, name), surface_space)
    canopy_height = ClimaLand.Canopy.clm_canopy_height(surface_space)
    tower = column_field(tower_height, surface_space)
    atmos_h = @. max(tower, canopy_height + FT(MIN_CANOPY_CLEARANCE))
    atmos = ClimaLand.PrescribedAtmosphere(
        column_input(:P_liq),
        column_input(:P_snow),
        column_input(:T),
        column_input(:u),
        column_input(:q),
        column_input(:P),
        SPINUP_START,
        atmos_h,
        toml_dict;
        c_co2 = column_input(:c_co2),
    )

    earth_param_set = LP.LandParameters(toml_dict)
    coordinates = ClimaCore.Fields.coordinate_field(surface_space)
    cos_zenith_angle =
        (t, s) -> ClimaLand.default_cos_zenith_angle(
            t,
            s;
            latitude = coordinates.lat,
            longitude = coordinates.long,
            insol_params = earth_param_set.insol_params,
        )
    radiation = ClimaLand.PrescribedRadiativeFluxes(
        FT,
        column_input(:SW_d),
        column_input(:LW_d),
        SPINUP_START;
        cosθs = cos_zenith_angle,
        toml_dict,
    )
    return (; atmos, radiation)
end

"""
    setup_simulation(toml_dict, sites, site_data; stop_date = STOP_DATE)

Return the simulation of the default land model on one column per site of
`sites` (see `read_sites`), forced by the site data loaded from
`SITE_DATA_FILE`, and the `DictWriter` collecting the monthly means of
`SHORT_NAMES`.
"""
function setup_simulation(toml_dict, sites, site_data; stop_date = STOP_DATE)
    start_date = SPINUP_START
    domain = ClimaLand.Domains.ColumnEnsemble(;
        zlim = (FT(-15), FT(0)),
        nelements = 15,
        dz_tuple = (FT(3), FT(0.05)),
        longlat = [FT.(ll) for ll in sites.longlat],
    )
    surface_space = domain.space.surface
    check_column_order(surface_space, sites.longlat)

    forcing = site_forcing(
        site_data.times,
        site_data.forcing,
        sites.tower_height,
        surface_space,
        toml_dict,
    )
    # The sites' own years precede MODIS and differ between sites
    LAI = ClimaLand.Canopy.prescribed_climatological_lai_modis(surface_space)
    model = ClimaLand.LandModel{FT}(
        forcing,
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
the forcing and the outputs are matched to the sites by position.
"""
function check_column_order(surface_space, longlat)
    coordinates = ClimaCore.Fields.coordinate_field(surface_space)
    longs = vec(Array(ClimaCore.Fields.field2array(coordinates.long)))
    lats = vec(Array(ClimaCore.Fields.field2array(coordinates.lat)))
    (
        longs ≈ [FT(lon) for (lon, _) in longlat] && lats ≈ [FT(lat) for (_, lat) in longlat]
    ) || error("The columns are not in the order of `longlat`")
    return nothing
end

"""
    calibration_months()

Return the first days of the months compared to observations.
"""
calibration_months() = collect(
    Dates.Date(CALIBRATION_START):Dates.Month(1):(Dates.Date(
        STOP_DATE,
    ) - Dates.Day(1)),
)

"""
    simulated_monthly_means(writer)

Return a dictionary mapping each of `SHORT_NAMES` to the `months × columns`
matrix of its monthly means in `calibration_months()`, in observation units.
"""
function simulated_monthly_means(writer)
    return Dict(
        map(SHORT_NAMES) do short_name
            saved = writer["$(short_name)_1M_average"]
            # A monthly mean is saved at the start of the following month
            month_of(t) = Dates.Date(date(t)) - Dates.Month(1)
            time_of_month = Dict(month_of(t) => t for t in keys(saved))
            values = reduce(
                vcat,
                [
                    permutedims(
                        vec(Array(parent(saved[time_of_month[month]]))),
                    ) for month in calibration_months()
                ],
            )
            short_name => get(SIM_TO_OBS_UNITS, short_name, identity).(values)
        end,
    )
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
    run_sites(interface::FluxnetSitesInterface, parameter_files, path)

Run the sites of `interface` with the default parameters overridden by the
calibrated parameters in the TOML files `parameter_files`, and save the
simulated monthly means to `path`.
"""
function run_sites(interface::FluxnetSitesInterface, parameter_files, path)
    override_files = model_parameter_file.(parameter_files)
    toml_dict = LP.create_toml_dict(FT; override_files)
    site_data = JLD2.load(site_data_path(interface))
    simulation, writer = setup_simulation(
        toml_dict,
        read_sites(sites_path(interface)),
        (; times = site_data["times"], forcing = site_data["forcing"]),
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
    JLD2.save_object(path, simulated_monthly_means(writer))
    return nothing
end

function ClimaCalibrate.forward_model(
    interface::FluxnetSitesInterface,
    iteration,
    member,
)
    (; output_dir) = interface
    member_path =
        ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, member)
    run_sites(
        interface,
        [ClimaCalibrate.parameter_path(output_dir, iteration, member)],
        joinpath(member_path, MONTHLY_MEANS_FILE),
    )
    return nothing
end

"""
    flatten(monthly, masks)

Stack the calibrated entries (`masks`) of the monthly means of each of
`SHORT_NAMES` into a vector ordered as the observations.
"""
flatten(monthly, masks) =
    reduce(vcat, [vec(monthly[name])[masks[name]] for name in SHORT_NAMES])

function ClimaCalibrate.observation_map(
    interface::FluxnetSitesInterface,
    iteration,
)
    (; output_dir) = interface
    masks = JLD2.load(joinpath(output_dir, OBSERVATIONS_FILE), "masks")
    ekp = ClimaCalibrate.load_ekp_struct(output_dir, iteration)
    n_observations = sum(sum, values(masks))
    g_sites = map(1:EKP.get_N_ens(ekp)) do member
        path = joinpath(
            ClimaCalibrate.path_to_ensemble_member(
                output_dir,
                iteration,
                member,
            ),
            MONTHLY_MEANS_FILE,
        )
        if isfile(path)
            flatten(JLD2.load_object(path), masks)
        else
            @error "Member $member has no output, filling it with NaNs"
            fill(NaN, n_observations)
        end
    end
    return reduce(hcat, g_sites)
end
