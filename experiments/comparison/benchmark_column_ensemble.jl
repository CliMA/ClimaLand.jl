# Benchmark the new space + ClimaUtilities forcing path against the old column, on the
# `run_fluxnet.jl` model (see `fluxnet_model.jl`):
#
#   old : `Column`         + FLUXNET CSV -> 0-D `TimeVaryingInput(times, values)`
#   new : `ColumnEnsemble` + FLUXNET-style NetCDF -> `DataSource` ->
#         `MultiColumnDataHandler` -> `TimeVaryingInput`
#
# One configuration per process, so compilation, GC state, and ClimaUtilities' global open-file
# cache never leak between the points of a sweep. Each run reports SYPD and appends a row to
# `output/benchmark_column_ensemble_<device>.txt`.
#
# Run: julia --project=.buildkite experiments/comparison/benchmark_column_ensemble.jl MODE N [NSTEPS] [SITE_ID]
#
#   MODE    `old` (requires N = 1) or `new`
#   N       number of columns
#   NSTEPS  timed steps, default 960 (five simulated days at dt = 450)
#   SITE_ID default US-MOz; must be a site with a CSV and a ClimaLand site file
#
# The device comes from ClimaComms, so one script covers both:
#
#   export CLIMACOMMS_DEVICE="CPU"   CLIMACOMMS_CONTEXT="SINGLETON"
#   export CLIMACOMMS_DEVICE="CUDA"  CLIMACOMMS_CONTEXT="SINGLETON"
#
# Both domain constructors default to `device = ClimaComms.device()`, and
# `ClimaComms.@elapsed device` synchronizes, so a GPU timing measures completed kernels rather
# than launches. Results land in output/benchmark_column_ensemble_<cpu|gpu>.txt -- one file per
# device, never mixed. On the cluster, drive this through sweep_column_ensemble.sh, which
# handles slurm and `module load climacommon`.
#
# A sweep, baseline first:
#
#   julia --project=.buildkite experiments/comparison/benchmark_column_ensemble.jl old 1
#   for N in 1 2 4 8 16 32; do
#     julia --project=.buildkite experiments/comparison/benchmark_column_ensemble.jl new $N
#   done
#
# What this does and does not measure:
# - `old` is N = 1 only. The 0-D `TimeVaryingInput` carries no spatial information, so the
#   baseline is a single point rather than a curve.
# - The N columns are all the same site, offset in longitude. ClimaUtilities matches a file to a
#   column by the file's scalar longitude/latitude and errors when two files share a location,
#   so the columns have to be distinguishable; the offset is 10x the matching tolerance and
#   physically negligible. A run over genuinely different sites would also carry per-column
#   `atmos_h`, UTC offsets, and site geometry, none of which this benchmark exercises.
# - MODIS LAI is read through ClimaUtilities in both modes, but its snapshots change monthly, so
#   no LAI read lands inside a timed window of a few days. This measures the forcing handler.
# - The forcing loading path IS inside the timed window: the handler reads and regrids a new
#   snapshot every time the pair of data dates bracketing `t` changes, which is every
#   `data_dt / dt` steps (4 at 30-minute forcing and dt = 450). The reported
#   `forcing_crossings` is how many times that happened, so a window that skipped the loading
#   path is visible rather than assumed.

# The buildkite environment sets this globally and it selects which CUDA.jl allocator is used;
# the scripts in experiments/benchmarks drop it before timing anything, so this one does too.
delete!(ENV, "JULIA_CUDA_MEMORY_POOL")

import ClimaComms
ClimaComms.@import_required_backends
using Dates
using DelimitedFiles # loads the FluxnetSimulations extension

using ClimaLand
using ClimaLand.Domains: Column, ColumnEnsemble
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: step!
import ClimaLand.FluxnetSimulations as FluxnetSimulations
import ClimaUtilities.FileReaders: close_all_ncfiles

include(joinpath(@__DIR__, "fluxnet_model.jl"))
include(joinpath(@__DIR__, "fluxnet_netcdf.jl"))

const FT = Float64

length(ARGS) >= 2 || error(
    "Usage: benchmark_column_ensemble.jl MODE N [NSTEPS] [SITE_ID], \
     where MODE is `old` or `new`",
)
const MODE = ARGS[1]
MODE in ("old", "new") || error("MODE must be `old` or `new`, got `$MODE`")
const NCOLUMNS = parse(Int, ARGS[2])
NCOLUMNS >= 1 || error("N must be at least 1, got $NCOLUMNS")
MODE == "old" &&
    NCOLUMNS != 1 &&
    error("`old` runs a single `Column`; pass N = 1 (got $NCOLUMNS)")
const NSTEPS = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 960
NSTEPS >= 1 || error("NSTEPS must be at least 1, got $NSTEPS")
const SITE_ID = length(ARGS) >= 4 ? ARGS[4] : "US-MOz"

context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.AbstractCPUDevice ? "cpu" : "gpu"

toml_dict = LP.create_toml_dict(FT)
earth_param_set = LP.LandParameters(toml_dict)
thermo_params = LP.thermodynamic_parameters(earth_param_set)

site = fluxnet_site_setup(FT, SITE_ID)
(; dt, lat, long, zmin, zmax, nelements, dz_tuple, time_offset, atmos_h) = site

# ---------------------------------------------------------------------------
# Forcing datasets, generated once and reused.
# ---------------------------------------------------------------------------
# Each column needs its own location: the handler matches files to columns by the file's scalar
# longitude/latitude within `atol = 1e-3` degrees and errors when two files match the same
# column. The offset depends only on the column index, never on N, so a smaller sweep point
# reuses the files a larger one generated.
const COLUMN_LON_OFFSET = FT(1e-2)
column_longlat(i) = (long + (i - 1) * COLUMN_LON_OFFSET, lat)

function generate_met_paths(ncolumns)
    dir = mkpath(joinpath(@__DIR__, "output", "bench_met"))
    return map(1:ncolumns) do i
        path = joinpath(dir, "$(SITE_ID)_Met_col$(lpad(i, 3, '0')).nc")
        if !isfile(path)
            (col_long, col_lat) = column_longlat(i)
            write_fluxnet_netcdf(
                path,
                SITE_ID,
                col_lat,
                col_long,
                time_offset,
                thermo_params,
                FT,
            )
            @info "Generated forcing NetCDF" path col_long col_lat
        end
        path
    end
end

met_paths = MODE == "new" ? generate_met_paths(NCOLUMNS) : String[]

# ---------------------------------------------------------------------------
# The two builds. Everything except the domain and the forcing comes from `site`.
# ---------------------------------------------------------------------------
function build_sim()
    if MODE == "old"
        domain = Column(;
            zlim = (zmin, zmax),
            nelements = nelements,
            dz_tuple = dz_tuple,
            longlat = (long, lat),
        )
        forcing = FluxnetSimulations.prescribed_forcing_fluxnet(
            SITE_ID,
            lat,
            long,
            time_offset,
            atmos_h,
            site.start_date,
            toml_dict,
            FT,
        )
    else
        domain = ColumnEnsemble(;
            zlim = (zmin, zmax),
            nelements = nelements,
            dz_tuple = dz_tuple,
            longlat = [column_longlat(i) for i in 1:NCOLUMNS],
        )
        # A scalar `atmos_h` and UTC offset apply to every column, which is what replicating one
        # site calls for; distinct sites would pass a surface Field and a vector instead.
        forcing = FluxnetSimulations.prescribed_forcing_netcdf(
            met_paths,
            domain.space.surface,
            atmos_h,
            site.start_date,
            toml_dict,
            FT;
            hour_offset_from_UTC = time_offset,
        )
    end
    # Neither constructor is given a device: both default to `ClimaComms.device()`, which reads
    # CLIMACOMMS_DEVICE. Check rather than assume, since a domain that silently landed on the CPU
    # would still produce clean-looking timings.
    domain_device = ClimaComms.device(domain)
    domain_device == device || error(
        "The domain was built on $domain_device but the benchmark is timing $device",
    )
    return build_fluxnet_sim(
        domain,
        forcing,
        site,
        toml_dict;
        user_callbacks = (),
    )
end

# `ClimaComms.@elapsed` runs its body inside a closure, so a binding assigned there does not
# escape; the result comes back through a `Ref`.
function timed(f, device)
    result = Ref{Any}(nothing)
    seconds = ClimaComms.@elapsed device result[] = f()
    return (result[], seconds)
end

# Stepping through a function keeps the simulation's type out of the global scope, so the timed
# loop does not pay a dynamic dispatch per step.
function step_n!(sim, n)
    for _ in 1:n
        step!(sim)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Warmup. The NetCDF path only touches disk when the two data dates bracketing `t` change, and
# the handler caches exactly those two snapshots, so the warmup has to cross several of those
# boundaries to compile the read, the composed Precip/VPD functions, and the cache eviction.
# ---------------------------------------------------------------------------
data_dt = fluxnet_data_dt(SITE_ID)
warmup_steps = max(4, ceil(Int, 3 * data_dt / dt))

@info "Warmup" mode = MODE ncolumns = NCOLUMNS site = SITE_ID data_dt warmup_steps dt device device_suffix
let (warmup, warmup_setup_s) = timed(build_sim, device)
    warmup_step_s = ClimaComms.@elapsed device step_n!(warmup, warmup_steps)
    @info "Warmup complete (includes compilation)" warmup_setup_s warmup_step_s
end
GC.gc()
# Drop ClimaUtilities' global open-file cache so the timed build reopens its datasets cold.
close_all_ncfiles()
GC.gc()

# ---------------------------------------------------------------------------
# The measured run.
# ---------------------------------------------------------------------------
(sim, setup_s) = timed(build_sim, device)
first_step_s = ClimaComms.@elapsed device step_n!(sim, 1)

# Timed in chunks: one aggregate number cannot show a GC pause or a warmup tail, and the spread
# across the window is what says whether the mean is worth quoting. SYPD uses the total.
const NCHUNKS = min(4, NSTEPS)
chunk_steps = fill(NSTEPS ÷ NCHUNKS, NCHUNKS)
chunk_steps[end] += NSTEPS - sum(chunk_steps)
chunk_s = map(n -> ClimaComms.@elapsed(device, step_n!(sim, n)), chunk_steps)
@assert sum(isnan.(sim._integrator.u)) == 0

wall_s = sum(chunk_s)
s_per_step = wall_s / NSTEPS
(s_per_step_min, s_per_step_max) = extrema(chunk_s ./ chunk_steps)

# Each crossing of the pair of data dates bracketing `t` makes every forcing variable read and
# regrid a new snapshot for every column, because the handler's LRU holds exactly those two
# brackets. This is therefore how many times the timed window exercised the loading path.
forcing_crossings = NSTEPS * dt / data_dt
forcing_crossings >= 2 || @warn "The timed window spans fewer than two forcing intervals, so \
                                 the data-loading path is barely exercised. Raise NSTEPS." forcing_crossings data_dt

# SYPD is simulated years of model time per wall-clock day. Every column advances on the same
# `dt`, so the model time covered does not depend on N; `column_sypd` is the throughput that
# does. 365.25 days/year matches `ClimaUtilities.OnlineLogging.sypd_str_from_ssps`.
sim_seconds = NSTEPS * dt
sypd = sim_seconds / wall_s / 365.25
column_sypd = sypd * NCOLUMNS

@info "Per-chunk wall time (s)" chunk_steps chunk_s
@info "Benchmark" mode = MODE ncolumns = NCOLUMNS nsteps = NSTEPS forcing_crossings dt setup_s first_step_s wall_s s_per_step s_per_step_min s_per_step_max sypd column_sypd

# ---------------------------------------------------------------------------
# Append one row per invocation, so a sweep accumulates in one file.
# ---------------------------------------------------------------------------
const RESULT_COLUMNS = (
    "date",
    "git_sha",
    "host",
    "device",
    "nthreads",
    "mode",
    "ncolumns",
    "site",
    "nelements",
    "dt",
    "nsteps",
    "warmup_steps",
    "forcing_crossings",
    "setup_s",
    "first_step_s",
    "wall_s",
    "s_per_step",
    "s_per_step_min",
    "s_per_step_max",
    "sypd",
    "column_sypd",
)

git_sha() =
    try
        readchomp(`git -C $(@__DIR__) rev-parse --short HEAD`)
    catch
        "unknown"
    end

row = (
    Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS"),
    git_sha(),
    Base.Libc.gethostname(),
    device_suffix,
    Threads.nthreads(),
    MODE,
    NCOLUMNS,
    SITE_ID,
    nelements,
    dt,
    NSTEPS,
    warmup_steps,
    round(forcing_crossings; sigdigits = 5),
    round(setup_s; sigdigits = 5),
    round(first_step_s; sigdigits = 5),
    round(wall_s; sigdigits = 5),
    round(s_per_step; sigdigits = 5),
    round(s_per_step_min; sigdigits = 5),
    round(s_per_step_max; sigdigits = 5),
    round(sypd; sigdigits = 5),
    round(column_sypd; sigdigits = 5),
)

results_path = joinpath(
    mkpath(joinpath(@__DIR__, "output")),
    "benchmark_column_ensemble_$(device_suffix).txt",
)
expected_header = join(RESULT_COLUMNS, "\t")
write_header = !isfile(results_path) || filesize(results_path) == 0
# Appending under a header from an older column set would silently misalign every later row.
write_header || open(readline, results_path) == expected_header || error(
    "$results_path was written with a different set of columns; move it aside and rerun.",
)
open(results_path, "a") do io
    write_header && println(io, expected_header)
    println(io, join(row, "\t"))
end
@info "Appended result" results_path

close_all_ncfiles()
