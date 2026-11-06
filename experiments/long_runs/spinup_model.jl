# # Global spin-up of the land model

# Spins up the soil, snow, canopy and lake state of the global land model by
# repeating the ERA5 forcing of 1979-1998. The simulation calendar starts on
# 2000-03-01 and each 20-year cycle sees the same forcing, so the drift between
# consecutive cycles measures how far the state is from equilibrium.

# The run starts from the latest checkpoint under `SPINUP_ROOT`; without one,
# from the checkpoint `RESTART_FILE`; without that, from the initial conditions
# file `IC_PATH` (by default the `rosetta_spunup_ic` artifact). Each invocation
# advances at most `CHUNK_YEARS` years (or `CHUNK_MONTHS` months), writing a
# checkpoint every `CHECKPOINT_MONTHS` months, until `N_CYCLES` cycles are done;
# it then writes `DONE_<year>` in `SPINUP_ROOT`.

# The monthly means of the last February (`swc`, `si`, `sie`, `swe`, `lwp`) are
# the input of `rosetta_spunup_ic/create_artifacts.jl` in ClimaArtifacts.

# Environment variables: `SPINUP_ROOT`, `N_CYCLES` (3), `CHUNK_YEARS` (20),
# `CHUNK_MONTHS`, `CHECKPOINT_MONTHS` (12), `IC_PATH`, `RESTART_FILE`,
# `PARAM_OVERRIDE_TOML`, `ALLOW_CPU`.

import ClimaComms
ClimaComms.@import_required_backends
import ClimaCore: Fields
import ClimaCore.InputOutput: HDF5
import ClimaUtilities
import ClimaUtilities.TimeManager: ITime, date
import ClimaParams as CP
using ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
using Dates

const FT = Float64
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device isa ClimaComms.CPUSingleThreaded &&
    !haskey(ENV, "ALLOW_CPU") &&
    error("The spin-up must run on a GPU; set CLIMACOMMS_DEVICE=CUDA.")

const START_DATE = DateTime(2000, 3, 1)
const CYCLE_YEARS = 20
const N_CYCLES = parse(Int, get(ENV, "N_CYCLES", "3"))
const STOP_DATE = START_DATE + Year(CYCLE_YEARS * N_CYCLES)
const CHUNK =
    haskey(ENV, "CHUNK_MONTHS") ? Month(parse(Int, ENV["CHUNK_MONTHS"])) :
    Year(parse(Int, get(ENV, "CHUNK_YEARS", "20")))
const CHECKPOINT_PERIOD = Month(parse(Int, get(ENV, "CHECKPOINT_MONTHS", "12")))
const Δt = 900.0

root_path = get(ENV, "SPINUP_ROOT", "snowy_land_spinup")
diagnostics_root = joinpath(root_path, "global_diagnostics")
mkpath(diagnostics_root)

function setup_model(::Type{FT}, Δt, domain, toml_dict) where {FT}
    surface_space = domain.space.surface
    # Periodic in time, so the 20 years repeat over the simulation
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        DateTime(1979),
        DateTime(1998),
        surface_space,
        toml_dict,
        FT;
        max_wind_speed = 25.0,
        context,
    )
    LAI = ClimaLand.Canopy.prescribed_climatological_lai_modis(surface_space)
    return LandModel{FT}(
        (; atmos, radiation),
        LAI,
        toml_dict,
        domain,
        Δt;
        prognostic_land_components = (:canopy, :lake, :snow, :soil),
    )
end

domain =
    ClimaLand.Domains.global_box_domain(FT; context, mask_threshold = FT(0.99))
override_files =
    haskey(ENV, "PARAM_OVERRIDE_TOML") ? [ENV["PARAM_OVERRIDE_TOML"]] : String[]
toml_dict = LP.create_toml_dict(FT; override_files)
model = setup_model(FT, Δt, domain, toml_dict)

# ClimaCore cannot rebuild the latitude-longitude box domain from a checkpoint
# ("Invalid coord type LatPoint"), so a restart copies the arrays saved in the
# checkpoint into the state of this model, which lives on the same spaces.
"""
    checkpoint_time(file)

Return the time of the checkpoint `file`, as an `ITime` with epoch `START_DATE`.
"""
function checkpoint_time(file)
    return HDF5.h5open(file, "r") do h5
        epoch = HDF5.read_attribute(h5, "epoch")
        epoch == string(START_DATE) ||
            error("$file has epoch $epoch, expected $START_DATE")
        seconds = round(Int, HDF5.read_attribute(h5, "time"))
        return ITime(seconds, Second(1), START_DATE)
    end
end

"""
    copy_checkpoint!(Y, file)

Copy the state saved in the checkpoint `file` into `Y`.
"""
copy_checkpoint!(Y, file) =
    HDF5.h5open(h5 -> copy_fields!(Y, h5["fields/Y"]), file, "r")

function copy_fields!(Y::Fields.FieldVector, group)
    for name in propertynames(Y)
        copy_fields!(getproperty(Y, name), group[string(name)])
    end
    return nothing
end

function copy_fields!(field::Fields.Field, dataset)
    data = read(dataset)
    size(data) == size(parent(field)) || error(
        "$(HDF5.name(dataset)) has size $(size(data)), expected $(size(parent(field)))",
    )
    copyto!(parent(field), data)
    return nothing
end

restart_file = ClimaLand.find_restart(diagnostics_root)
if isnothing(restart_file) && haskey(ENV, "RESTART_FILE")
    restart_file = ENV["RESTART_FILE"]
end
if isnothing(restart_file)
    # Never fall back to the initial conditions in the middle of a spin-up
    isempty(filter(startswith("output_"), readdir(diagnostics_root))) || error(
        "No checkpoint in $diagnostics_root, which has output folders: remove them to start from the initial conditions",
    )
    ic_path = get(
        ENV,
        "IC_PATH",
        ClimaLand.Artifacts.rosetta_spunup_ic_path(; context),
    )
    @info "Starting from the initial conditions in $ic_path"
    t_start = ITime(0, Second(1), START_DATE)
    set_ic! =
        ClimaLand.Simulations.make_set_initial_state_from_file(ic_path, model)
else
    t_start = checkpoint_time(restart_file)
    @info "Restarting from $restart_file ($(date(t_start)))"
    set_ic! = (Y, p, t, model) -> copy_checkpoint!(Y, restart_file)
end

chunk_start = date(t_start)
chunk_stop = min(STOP_DATE, chunk_start + CHUNK)
done_file = joinpath(root_path, "DONE_$(year(STOP_DATE))")
if chunk_start >= STOP_DATE
    @info "Already at $chunk_start, nothing to do (stop date $STOP_DATE)"
    touch(done_file)
    exit(0)
end

outdir =
    ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_root)
t0, tf, Δt_itime = promote(
    t_start,
    ITime(Dates.value(Second(chunk_stop - START_DATE)); epoch = START_DATE),
    ITime(Δt; epoch = START_DATE),
)

diagnostics = ClimaLand.default_diagnostics(
    model,
    START_DATE,
    outdir;
    reduction_period = :monthly,
    output_vars = [
        "tsoil",
        "sie",
        "swc",
        "si",
        "lwp",
        "tair",
        "swe",
        "lhf",
        "shf",
    ],
    date_last = chunk_start,
)

# Abort on NaN rather than warn, so that a failed chunk stops the job chain
mask = ClimaLand.Domains.landsea_mask(domain)
nan_count(field) = mapreduce(
    (s, m) -> m != 0 && isnan(s),
    Base.add_sum,
    parent(field),
    parent(mask),
)
top(field) = ClimaLand.Domains.top_center_to_surface(field)
function abort_on_nan(integrator)
    Y = integrator.u
    n =
        nan_count(top(Y.soil.ϑ_l)) +
        nan_count(top(Y.soil.ρe_int)) +
        nan_count(Y.snow.U) +
        nan_count(Y.canopy.energy.T)
    n > 0 && error("$n NaNs in the land state on $(date(integrator.t))")
    return nothing
end

simulation = LandSimulation(
    t0,
    tf,
    Δt_itime,
    model;
    outdir,
    set_ic!,
    user_callbacks = (
        ClimaLand.IntervalBasedCallback(Month(1), t0, Δt_itime, abort_on_nan),
        ClimaLand.ReportCallback(div((tf - t0), 10), t0),
        ClimaLand.CheckpointCallback(
            CHECKPOINT_PERIOD,
            outdir,
            t0;
            model,
            dt = Δt_itime,
        ),
    ),
    diagnostics,
)
@info "Spin-up of the global land model"
@info "Device: $device"
@info "Resolution: $(domain.nelements)"
@info "Chunk: $chunk_start to $chunk_stop (stop date $STOP_DATE)"
CP.log_parameter_information(toml_dict, joinpath(root_path, "parameters.toml"))
solve!(simulation)

# A checkpoint at the end of the chunk, in case the callback did not fire there
ClimaLand.save_checkpoint(
    simulation._integrator.u,
    simulation._integrator.t,
    outdir;
    model,
)
date(simulation._integrator.t) >= STOP_DATE && touch(done_file)
@info "Chunk done at $(date(simulation._integrator.t))"
