# Plot the column-ensemble benchmark sweeps written by `benchmark_column_ensemble.jl`.
#
# Reads every per-device results file, keeps the most recent row per (device, mode, N), prints the
# table it used, and writes a three-panel PNG:
#
#   1. rate       -- SYPD, simulated years of model time per wall day. Falls with N, because all
#                    N columns advance on one shared `dt`.
#   2. throughput -- SYPD x N, column-years per wall day. This is the SYPD normalized by the
#                    column count: flat means N columns cost exactly N times one column, rising
#                    means the ensemble amortizes per-step overhead.
#   3. build cost -- setup seconds, where opening and matching N forcing files shows up.
#
# CPU and GPU are separate files and are drawn together: colour is the device, and the line style
# is the configuration -- solid for the `ColumnEnsemble`, dashed for the old single `Column`
# reference. That reference is read as N of those runs with one per processor, so it is flat in
# the rate and build panels (what one processor does) and grows as N x that rate in the throughput
# panel (ideal embarrassingly-parallel scaling). Note the asymmetry that makes it a demanding
# reference -- the reference gets N processors, the ensemble curve gets one.
#
# Run: julia --project=.buildkite experiments/comparison/plot_column_ensemble_benchmark.jl [SITE_ID] [RESULTS_FILE...]
#
# With no RESULTS_FILE, both output/benchmark_column_ensemble_cpu.txt and _gpu.txt are used when
# they exist, so the same command works before and after a GPU sweep.

using CairoMakie
using DelimitedFiles
using Printf

const SITE_ID = length(ARGS) >= 1 ? ARGS[1] : "US-MOz"
const RESULTS_PATHS =
    length(ARGS) >= 2 ? ARGS[2:end] :
    filter(
        isfile,
        [
            joinpath(
                @__DIR__,
                "output",
                "benchmark_column_ensemble_$(device).txt",
            ) for device in ("cpu", "gpu")
        ],
    )

isempty(RESULTS_PATHS) && error(
    "No results files found under $(joinpath(@__DIR__, "output")); run \
     sweep_column_ensemble.sh first, or pass paths as arguments",
)

# Chart chrome and the categorical slots, from the shared palette.
const SURFACE = "#fcfcfb"
const INK = "#0b0b0b"
const SECONDARY = "#52514e"
const MUTED = "#898781"
const GRID = "#e1e0d9"
const AXIS = "#c3c2b7"
const DEVICE_COLORS = Dict("cpu" => "#2a78d6", "gpu" => "#eb6834")
const FALLBACK_COLOR = "#1baf7a"
device_color(device) = get(DEVICE_COLORS, device, FALLBACK_COLOR)

# ---------------------------------------------------------------------------
# Read the results. Each file carries its own header, so a file written before a column was
# added still loads.
# ---------------------------------------------------------------------------
function read_results(path)
    (rows, header) = readdlm(path, '\t'; header = true)
    index = Dict(String(name) => i for (i, name) in enumerate(vec(header)))
    cell(row, name) = rows[row, index[name]]
    as_int(x) = x isa AbstractString ? parse(Int, x) : round(Int, x)
    as_float(x) = x isa AbstractString ? parse(Float64, x) : Float64(x)
    float_or_nan(row, name) =
        haskey(index, name) ? as_float(cell(row, name)) : NaN
    return [
        (;
            device = String(cell(row, "device")),
            mode = String(cell(row, "mode")),
            site = String(cell(row, "site")),
            ncolumns = as_int(cell(row, "ncolumns")),
            nsteps = as_int(cell(row, "nsteps")),
            wall_s = as_float(cell(row, "wall_s")),
            s_per_step = as_float(cell(row, "s_per_step")),
            s_per_step_min = float_or_nan(row, "s_per_step_min"),
            s_per_step_max = float_or_nan(row, "s_per_step_max"),
            forcing_crossings = float_or_nan(row, "forcing_crossings"),
            setup_s = as_float(cell(row, "setup_s")),
            sypd = as_float(cell(row, "sypd")),
            column_sypd = as_float(cell(row, "column_sypd")),
        ) for row in axes(rows, 1)
    ]
end

# A sweep can be re-run, so the files hold repeats; the latest row for a (device, mode, N) wins.
records = Dict{Tuple{String, String, Int}, NamedTuple}()
superseded = 0
for path in RESULTS_PATHS, point in read_results(path)
    point.site == SITE_ID || continue
    key = (point.device, point.mode, point.ncolumns)
    haskey(records, key) && (superseded += 1)
    records[key] = point
end

isempty(records) && error(
    "No rows for site $SITE_ID in $(join(RESULTS_PATHS, ", "))",
)
by_columns = points -> sort(points; by = p -> p.ncolumns)
devices = sort(unique(p.device for p in values(records)))
series_for(device, mode) =
    by_columns([
        p for p in values(records) if p.device == device && p.mode == mode
    ])
all_points = sort(
    collect(values(records));
    by = p -> (p.device, p.mode, p.ncolumns),
)

# ---------------------------------------------------------------------------
# The numbers behind the figure.
# ---------------------------------------------------------------------------
println()
println("$(join(RESULTS_PATHS, "\n"))")
println("site $SITE_ID, devices: $(join(devices, ", "))")
superseded > 0 && println(
    "$superseded earlier row(s) superseded by a later run of the same (device, mode, N)",
)
@printf(
    "%-4s %-5s %4s %7s %9s %9s %11s %11s %11s %9s %11s %8s\n",
    "dev",
    "mode",
    "N",
    "nsteps",
    "crossings",
    "wall_s",
    "s_per_step",
    "s/step_min",
    "s/step_max",
    "sypd",
    "sypd*N",
    "setup_s"
)
for p in all_points
    @printf(
        "%-4s %-5s %4d %7d %9.1f %9.2f %11.5f %11.5f %11.5f %9.4f %11.4f %8.2f\n",
        p.device,
        p.mode,
        p.ncolumns,
        p.nsteps,
        p.forcing_crossings,
        p.wall_s,
        p.s_per_step,
        p.s_per_step_min,
        p.s_per_step_max,
        p.sypd,
        p.column_sypd,
        p.setup_s
    )
end
println()

# ---------------------------------------------------------------------------
# The figure.
# ---------------------------------------------------------------------------
ns = sort(unique(p.ncolumns for p in all_points))
xscale = length(ns) > 1 ? log2 : identity

function panel(fig, col; title, ylabel, yscale = identity)
    return Axis(
        fig[1, col];
        title,
        titlealign = :left,
        titlecolor = INK,
        titlesize = 14,
        titlefont = :bold,
        xlabel = "columns (N)",
        ylabel,
        xlabelcolor = SECONDARY,
        ylabelcolor = SECONDARY,
        xlabelsize = 12,
        ylabelsize = 12,
        xticks = (ns, string.(ns)),
        xscale,
        yscale,
        xticklabelcolor = MUTED,
        yticklabelcolor = MUTED,
        xticklabelsize = 11,
        yticklabelsize = 11,
        xtickcolor = AXIS,
        ytickcolor = AXIS,
        xgridcolor = GRID,
        ygridcolor = GRID,
        xgridwidth = 1,
        ygridwidth = 1,
        leftspinecolor = AXIS,
        bottomspinecolor = AXIS,
        topspinevisible = false,
        rightspinevisible = false,
        backgroundcolor = SURFACE,
    )
end

label_value(x) =
    abs(x) >= 100 ? @sprintf("%.0f", x) :
    abs(x) >= 10 ? @sprintf("%.1f", x) : @sprintf("%.3g", x)

# Only the endpoint is labelled: a number on every point is noise.
function annotate_end!(ax, x, y; align, offset)
    text!(
        ax,
        x,
        y;
        text = label_value(y),
        color = SECONDARY,
        fontsize = 11,
        align,
        offset,
    )
    return nothing
end

function draw_series!(ax, points, value; color, label)
    isempty(points) && return nothing
    xs = [p.ncolumns for p in points]
    ys = [value(p) for p in points]
    if length(xs) > 1
        lines!(ax, xs, ys; color, linewidth = 2, label)
        scatter!(
            ax,
            xs,
            ys;
            color,
            markersize = 10,
            strokecolor = SURFACE,
            strokewidth = 2,
        )
    else
        scatter!(
            ax,
            xs,
            ys;
            color,
            markersize = 10,
            strokecolor = SURFACE,
            strokewidth = 2,
            label,
        )
    end
    annotate_end!(
        ax,
        xs[end],
        ys[end];
        align = (:center, :bottom),
        offset = (0, 9),
    )
    return nothing
end

# The old column is a single measurement, but it reads as the level to beat, so it spans the
# panel. `scaled` gives the throughput reading: N of those runs, one per processor, deliver N
# times one run's column-years per wall day. Unscaled is the per-processor reading, which does
# not depend on how many processors are running.
function draw_reference!(ax, points, value; scaled, color, label)
    isempty(points) && return nothing
    level = value(first(points))
    if scaled
        ys = level .* ns
        lines!(ax, ns, ys; color, linewidth = 2, linestyle = :dash, label)
        annotate_end!(
            ax,
            ns[end],
            ys[end];
            align = (:right, :top),
            offset = (-4, -6),
        )
    else
        hlines!(ax, [level]; color, linewidth = 2, linestyle = :dash, label)
    end
    return nothing
end

fig = Figure(;
    size = (1240, 460),
    backgroundcolor = SURFACE,
    figure_padding = 26,
)

steps = join(sort(unique(p.nsteps for p in all_points)), ", ")
Label(
    fig[0, 1:3],
    "Column ensemble vs single column — $SITE_ID, $steps timed steps";
    color = INK,
    fontsize = 17,
    font = :bold,
    halign = :left,
    padding = (0, 0, 4, 0),
)

# The reference scales with N only where the measure is a throughput; the log y on that panel
# keeps both it and the ensemble curve readable when they diverge by more than 10x.
panels = (
    (
        panel(fig, 1; title = "Rate", ylabel = "SYPD", yscale = log10),
        p -> p.sypd,
        false,
    ),
    (
        panel(
            fig,
            2;
            title = "Throughput — column-years / wall day",
            ylabel = "SYPD × N",
            yscale = log10,
        ),
        p -> p.column_sypd,
        true,
    ),
    (
        panel(fig, 3; title = "Build cost", ylabel = "setup (s)"),
        p -> p.setup_s,
        false,
    ),
)

for (ax, value, scaled) in panels, device in devices
    color = device_color(device)
    draw_series!(
        ax,
        series_for(device, "new"),
        value;
        color,
        label = "$(uppercase(device)) — ColumnEnsemble",
    )
    draw_reference!(
        ax,
        series_for(device, "old"),
        value;
        scaled,
        color,
        label = "$(uppercase(device)) — Column, one per processor",
    )
end

Legend(
    fig[2, 1:3],
    first(first(panels));
    orientation = :horizontal,
    framevisible = false,
    labelcolor = SECONDARY,
    labelsize = 12,
    nbanks = length(devices) > 1 ? 2 : 1,
)

out_path = joinpath(
    @__DIR__,
    "output",
    "benchmark_column_ensemble_$(SITE_ID).png",
)
save(out_path, fig)
println("Wrote $out_path")
