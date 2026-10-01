# Energy and carbon RMSE boxplots: compare ClimaLand against an ILAMB land-hist
# cohort. "Other-model" RMSE values are inlined from the global table of the
# ILAMB land-hist dashboards (https://www.ilamb.org/land-hist/): CLM,
# ISBA-CTRIP and JSBACH (no LWup), each forced by CRUJRA, GSWP3 and Princeton,
# in that order; the multi-model "Mean-*" rows are left out. ClimaLand RMSE is
# computed at runtime from the diagnostics directory against the same
# benchmark as the cohort, with ILAMB's default RMSE definition.

# Energy benchmarks are read from the ILAMB `DATA` tree, as
# `<var>/<source>/<file>.nc`, rooted at `ENV["ILAMB_ROOT"]` (the convention of
# the ILAMB package) or at the copy on the CliMA cluster.
const _DEFAULT_ILAMB_ROOT = "/net/sampo/data1/ilamb"

# Panel definitions: sim short_name, benchmark name, benchmark file, factor
# converting the benchmark to W m⁻², and cohort RMSEs.
const _ENERGY_PANELS = (
    (
        title = "H",
        sim_short_name = "shf",
        bench = "FLUXCOM",
        obs_path = ("hfss", "FLUXCOM", "hfss.nc"),
        obs_to_sim_units = 1e6 / 86400, # MJ m⁻² day⁻¹ to W m⁻²
        others = [18.7, 16.3, 16.3, 20.8, 21.0, 18.9, 27.3, 27.3, 27.4],
    ),
    (
        title = "LE",
        sim_short_name = "lhf",
        bench = "FLUXCOM",
        obs_path = ("hfls", "FLUXCOM", "hfls.nc"),
        obs_to_sim_units = 1e6 / 86400,
        others = [16.4, 16.8, 19.3, 17.9, 18.9, 19.4, 25.8, 24.3, 24.2],
    ),
    (
        title = "SWup",
        sim_short_name = "swu",
        bench = "CERES",
        obs_path = ("rsus", "CERESed4.2", "rsus.nc"),
        obs_to_sim_units = 1.0,
        others = [11.1, 10.5, 11.0, 12.7, 11.7, 12.1, 12.9, 11.5, 12.3],
    ),
    (
        title = "LWup",
        sim_short_name = "lwu",
        bench = "CERES",
        obs_path = ("rlus", "CERESed4.2", "rlus.nc"),
        obs_to_sim_units = 1.0,
        others = [14.1, 14.7, 14.4, 14.5, 13.0, 13.2],
    ),
)

const _CARBON_PANELS = (
    (
        title = "GPP",
        sim_short_name = "gpp",
        bench = "FLUXCOM",
        others = [1.85, 1.73, 2.08, 1.5, 1.49, 1.81, 2.06, 1.85, 1.97],
    ),
    (
        title = "ER",
        sim_short_name = "er",
        bench = "FLUXCOM",
        others = [1.57, 1.32, 1.89, 1.43, 1.23, 1.75, 1.81, 1.53, 1.74],
    ),
)

"""
    _monthly_climatology_maps(var)

Return a 12-element vector of lon-lat arrays, the mean of `var` over all
samples of each calendar month. `NaN` samples are skipped; cells with no valid
sample in a month are `NaN`.
"""
function _monthly_climatology_maps(var)
    months = Dates.month.(ClimaAnalysis.dates(var))
    grid_size = size(_mask_template(var).data)
    return map(1:12) do m
        total = zeros(grid_size)
        count = zeros(grid_size)
        for t in ClimaAnalysis.times(var)[months .== m]
            x = ClimaAnalysis.slice(var, time = t).data
            valid = isfinite.(x)
            total[valid] .+= x[valid]
            count .+= valid
        end
        total ./ count
    end
end

"""
    _ilamb_cycle_rmse(sim_var, obs_var, mask_fn)

Return ILAMB's default global RMSE (`rmse_score_basis = "cycle"`) of `sim_var`
against `obs_var`, two monthly `OutputVar`s on the same grid and times. At each
cell, the RMSE is taken over the 12 months of the mean annual cycle of each
field; the global value is the area-weighted mean of that map over the cells
kept by `mask_fn`. Unlike the RMSE of time means, this includes seasonal-cycle
errors.
"""
function _ilamb_cycle_rmse(sim_var, obs_var, mask_fn)
    template = _mask_template(sim_var)
    sum_sq = zeros(size(template.data))
    n_months = zeros(size(template.data))
    for (sim_m, obs_m) in zip(
        _monthly_climatology_maps(sim_var),
        _monthly_climatology_maps(obs_var),
    )
        sq = (sim_m .- obs_m) .^ 2
        valid = isfinite.(sq)
        sum_sq[valid] .+= sq[valid]
        n_months .+= valid
    end
    rmse_map = ClimaAnalysis.remake(template; data = sqrt.(sum_sq ./ n_months))
    return ClimaAnalysis.weighted_average_lonlat(mask_fn(rmse_map)).data[]
end

"""
    _ilamb_energy_benchmark(panel)

Load the benchmark of energy `panel` from the ILAMB `DATA` tree as an
`OutputVar` in W m⁻² that follows CliMA conventions. Return `nothing` if the
file is not found.
"""
function _ilamb_energy_benchmark(panel)
    ilamb_root = get(ENV, "ILAMB_ROOT", _DEFAULT_ILAMB_ROOT)
    path = joinpath(ilamb_root, "DATA", panel.obs_path...)
    if !isfile(path)
        @warn "ILAMB benchmark $path not found; set ILAMB_ROOT to the ILAMB data root"
        return nothing
    end
    var = ClimaAnalysis.OutputVar(path, first(panel.obs_path))
    for (dim, dim_units) in (("lon", "degrees_east"), ("lat", "degrees_north"))
        ClimaAnalysis.dim_units(var, dim) == "degree" &&
            ClimaAnalysis.set_dim_units!(var, dim, dim_units)
    end
    ClimaAnalysis.transform_dates!(var, Dates.firstdayofmonth)
    replace!(var, missing => NaN)
    var = ClimaAnalysis.convert_units(
        var,
        "W m^-2",
        conversion_function = x -> x * panel.obs_to_sim_units,
    )
    ClimaAnalysis.set_short_name!(var, panel.sim_short_name)
    return _preprocess_var(var)
end

"""
    _boxplot_benchmark(panel)

Return the benchmark `OutputVar` of `panel` and a function of
`(sim_var, obs_var)` returning its mask function, or `nothing` if the benchmark
is unavailable. Energy benchmarks are masked to land cells where the benchmark
is defined; carbon benchmarks use the `ILAMBDataLoader` masks.
"""
function _boxplot_benchmark(panel)
    if hasproperty(panel, :obs_path)
        obs_var = _ilamb_energy_benchmark(panel)
        isnothing(obs_var) && return nothing
        make_mask =
            (sim_var, obs_var) -> begin
                valid_obs = ClimaAnalysis.make_lonlat_mask(
                    ClimaAnalysis.slice(
                        obs_var,
                        time = ClimaAnalysis.times(obs_var) |> first,
                    );
                    set_to_val = isnan,
                )
                return var -> ClimaAnalysis.apply_oceanmask(valid_obs(var))
            end
        return obs_var, make_mask
    end
    data_loader = ILAMBDataLoader()
    panel.sim_short_name in available_vars(data_loader) || return nothing
    return get(data_loader, panel.sim_short_name),
    get_mask_dict(data_loader)[panel.sim_short_name]
end

"""
    _boxplot_rmse(diagnostics_folder_path, panel; spin_up_months = 12)

Compute the global RMSE of ClimaLand `panel.sim_short_name` against the
benchmark of `panel`, the one the cohort values are scored against, as ILAMB
does (see `_ilamb_cycle_rmse`). The leaderboard figures keep the RMSE of time
means. Spinup is removed and both fields are windowed to their overlap and
resampled onto the simulation grid as in `compute_seasonal_leaderboard`.

Returns `NaN` if the variable is not in the simulation directory or the
benchmark is unavailable.
"""
function _boxplot_rmse(diagnostics_folder_path, panel; spin_up_months = 12)
    short_name = panel.sim_short_name
    sim_dir = ClimaAnalysis.SimDir(diagnostics_folder_path)
    short_name in ClimaAnalysis.available_vars(sim_dir) || return NaN
    benchmark = _boxplot_benchmark(panel)
    isnothing(benchmark) && return NaN
    obs_var, make_mask = benchmark

    sim_var = get(sim_dir, short_name)

    # Preprocess sim var to match conventions of data loaders
    sim_var = preprocess_sim_var(sim_var)

    ClimaAnalysis.set_reference_date!(obs_var, sim_var.attributes["start_date"])

    spinup_cutoff = spin_up_months * 31 * 86400.0
    if ClimaAnalysis.times(sim_var)[end] >= spinup_cutoff
        sim_var = ClimaAnalysis.window(sim_var, "time", left = spinup_cutoff)
    end

    sim_times = ClimaAnalysis.times(sim_var)
    obs_times = ClimaAnalysis.times(obs_var)
    min_time = maximum(first.((sim_times, obs_times)))
    max_time = minimum(last.((sim_times, obs_times)))
    sim_var =
        ClimaAnalysis.window(sim_var, "time", left = min_time, right = max_time)
    obs_var =
        ClimaAnalysis.window(obs_var, "time", left = min_time, right = max_time)

    obs_var = ClimaAnalysis.shift_longitude(obs_var, -180.0, 180.0)
    obs_var = ClimaAnalysis.resampled_as(obs_var, sim_var)

    return _ilamb_cycle_rmse(sim_var, obs_var, make_mask(sim_var, obs_var))
end

# Box-and-whisker statistics with Tukey-style 1.5*IQR fences clipped to data
# extrema.
function _bstats(vals)
    s = sort(collect(filter(isfinite, vals)))
    length(s) < 2 && return nothing
    q1 = Statistics.quantile(s, 0.25)
    med = Statistics.quantile(s, 0.50)
    q3 = Statistics.quantile(s, 0.75)
    iqr = q3 - q1
    lo = max(minimum(s), q1 - 1.5 * iqr)
    hi = min(maximum(s), q3 + 1.5 * iqr)
    return (; q1, med, q3, lo, hi)
end

# Draw one panel: cohort box + jittered model dots + ClimaLand current/prev dots.
function _draw_boxplot_panel!(
    fig,
    col,
    panel,
    rmse_current,
    rmse_prev,
    y_max,
    ylabel,
)
    xlabel = "vs $(panel.bench)\n(n=$(length(panel.others)) other models)"
    ax = CairoMakie.Axis(
        fig[1, col],
        title = panel.title,
        titlesize = 18,
        ylabel = ylabel,
        xlabel = xlabel,
        xlabelsize = 12,
        xlabelcolor = :gray,
        xticks = (Float64[], String[]),
        xgridvisible = false,
        limits = (0, 1, 0, y_max),
    )
    CairoMakie.hidexdecorations!(
        ax;
        ticks = true,
        ticklabels = true,
        label = false,
        grid = true,
    )

    cx = 0.40
    box_half = 0.085

    s = _bstats(panel.others)
    if !isnothing(s)
        # IQR box
        CairoMakie.poly!(
            ax,
            CairoMakie.Rect(cx - box_half, s.q1, 2 * box_half, s.q3 - s.q1);
            color = (:dodgerblue, 0.10),
            strokecolor = :dodgerblue,
            strokewidth = 1.5,
        )
        # Median
        CairoMakie.lines!(
            ax,
            [cx - box_half, cx + box_half],
            [s.med, s.med];
            color = :darkorange,
            linewidth = 3,
        )
        # Whiskers
        CairoMakie.lines!(
            ax,
            [cx, cx],
            [s.lo, s.q1];
            color = :dodgerblue,
            linewidth = 1.5,
        )
        CairoMakie.lines!(
            ax,
            [cx, cx],
            [s.q3, s.hi];
            color = :dodgerblue,
            linewidth = 1.5,
        )
        CairoMakie.lines!(
            ax,
            [cx - 0.04, cx + 0.04],
            [s.lo, s.lo];
            color = :dodgerblue,
            linewidth = 1.5,
        )
        CairoMakie.lines!(
            ax,
            [cx - 0.04, cx + 0.04],
            [s.hi, s.hi];
            color = :dodgerblue,
            linewidth = 1.5,
        )
    end

    # Jittered model dots (deterministic golden-ratio jitter for stable layout)
    for (i, v) in enumerate(panel.others)
        jx = cx + ((i * 0.618033) % 1) * 0.16 - 0.08
        CairoMakie.scatter!(
            ax,
            [jx],
            [v];
            color = (:dodgerblue, 0.4),
            strokecolor = :white,
            strokewidth = 0.7,
            markersize = 11,
        )
    end

    # ClimaLand dots offset to the right of the box
    dx = cx + box_half + 0.18

    if isfinite(rmse_prev) && isfinite(rmse_current)
        CairoMakie.lines!(
            ax,
            [dx, dx],
            [rmse_prev, rmse_current];
            color = :gray,
            linewidth = 1.5,
        )
    end
    if isfinite(rmse_prev)
        CairoMakie.scatter!(
            ax,
            [dx],
            [rmse_prev];
            color = :gray,
            strokecolor = :white,
            strokewidth = 2,
            markersize = 16,
        )
    end
    if isfinite(rmse_current)
        CairoMakie.scatter!(
            ax,
            [dx],
            [rmse_current];
            color = :firebrick,
            strokecolor = :white,
            strokewidth = 2,
            markersize = 16,
        )
    end
    if isfinite(rmse_prev) && isfinite(rmse_current) && rmse_prev > 0
        pct = round(Int, (rmse_prev - rmse_current) / rmse_prev * 100)
        sign = pct >= 0 ? "−" : "+"
        CairoMakie.text!(
            ax,
            dx + 0.04,
            (rmse_prev + rmse_current) / 2;
            text = "$(sign)$(abs(pct))%",
            color = pct >= 0 ? :forestgreen : :firebrick,
            fontsize = 14,
            font = :bold,
            align = (:left, :center),
        )
    end

    return ax
end

"""
    compute_rmse_boxplots(leaderboard_base_path,
                          diagnostics_folder_path;
                          prev_diagnostics_folder_path = nothing)

Generate `boxplot_rmse.png` in `leaderboard_base_path`: a single figure with
four energy panels (H, LE, SWup, LWup) and two carbon panels (GPP, ER) side
by side, separated by a slightly wider column gap because the two groups use
different units.

Each panel shows a boxplot of "other-model" RMSE values from the ILAMB
land-hist dashboards alongside the ClimaLand RMSE computed from the simulation
in `diagnostics_folder_path` (red dot). If `prev_diagnostics_folder_path` is
given, the equivalent RMSE from that earlier run is plotted as a gray dot,
with a percent-change label drawn between the two.

ClimaLand RMSE follows ILAMB's default definition (see `_ilamb_cycle_rmse`) and
is computed against the benchmark of the cohort: FLUXCOM for H and LE, CERES
EBAF for SWup and LWup (read from the ILAMB `DATA` tree, see
`_ilamb_energy_benchmark`), and ILAMB FLUXCOM for GPP and ER. Energy panels are
RMSE in W m⁻²; carbon panels are RMSE in g m⁻² day⁻¹.
"""
function compute_rmse_boxplots(
    leaderboard_base_path,
    diagnostics_folder_path;
    prev_diagnostics_folder_path = nothing,
)
    @info "Computing ClimaLand RMSE for energy/carbon boxplots"

    energy_panels = collect(_ENERGY_PANELS)
    carbon_panels = collect(_CARBON_PANELS)
    all_panels = vcat(energy_panels, carbon_panels)
    n_energy = length(energy_panels)

    rmse_current = Dict{String, Float64}()
    rmse_prev = Dict{String, Float64}()
    for p in all_panels
        rmse_current[p.sim_short_name] =
            _boxplot_rmse(diagnostics_folder_path, p)
        rmse_prev[p.sim_short_name] =
            isnothing(prev_diagnostics_folder_path) ? NaN :
            _boxplot_rmse(prev_diagnostics_folder_path, p)
    end

    function _group_y_max(panels)
        vals = Float64[]
        for p in panels
            append!(vals, p.others)
            push!(vals, rmse_current[p.sim_short_name])
            push!(vals, rmse_prev[p.sim_short_name])
        end
        finite = filter(isfinite, vals)
        return isempty(finite) ? 1.0 : maximum(finite) * 1.18
    end
    y_max_energy = _group_y_max(energy_panels)
    y_max_carbon = _group_y_max(carbon_panels)

    fig = CairoMakie.Figure(size = (260 * length(all_panels) + 280, 560))
    for (col, p) in enumerate(all_panels)
        is_energy = col <= n_energy
        y_max = is_energy ? y_max_energy : y_max_carbon
        ylabel = if col == 1
            "RMSE [W m⁻²]"
        elseif col == n_energy + 1
            "RMSE [g m⁻² day⁻¹]"
        else
            ""
        end
        _draw_boxplot_panel!(
            fig,
            col,
            p,
            rmse_current[p.sim_short_name],
            rmse_prev[p.sim_short_name],
            y_max,
            ylabel,
        )
    end
    # Widen the gap between the last energy panel and the first carbon panel
    # so the unit change reads as a deliberate split rather than another panel.
    CairoMakie.colgap!(fig.layout, n_energy, 50)

    CairoMakie.Label(
        fig[0, :],
        "ClimaLand RMSE — Global energy and carbon fluxes";
        fontsize = 18,
        font = :bold,
    )

    legend_elems = [
        CairoMakie.MarkerElement(
            color = (:dodgerblue, 0.4),
            marker = :circle,
            markersize = 11,
        ),
        CairoMakie.LineElement(color = :darkorange, linewidth = 3),
        CairoMakie.PolyElement(
            color = (:dodgerblue, 0.10),
            strokecolor = :dodgerblue,
            strokewidth = 1.5,
        ),
        CairoMakie.MarkerElement(
            color = :firebrick,
            marker = :circle,
            markersize = 14,
        ),
    ]
    legend_labels = [
        "Other models (ILAMB land-hist)",
        "Ensemble median",
        "IQR box",
        "ClimaLand · current",
    ]
    if !isnothing(prev_diagnostics_folder_path)
        push!(
            legend_elems,
            CairoMakie.MarkerElement(
                color = :gray,
                marker = :circle,
                markersize = 14,
            ),
        )
        push!(legend_labels, "ClimaLand · previous")
    end
    CairoMakie.Legend(
        fig[2, :],
        legend_elems,
        legend_labels;
        orientation = :horizontal,
        tellheight = true,
        framevisible = false,
        labelsize = 12,
    )

    disclaimer =
        "Disclaimer: For internal use only. A proper Model Intercomparison " *
        "Project (MIP) is beyond scope here, so this comparison is not " *
        "strictly fair: forcings differ across " *
        "models, tuning targets differ, and some variables are prescribed " *
        "rather than predicted (e.g., LAI in ClimaLand). Only a MIP would " *
        "enable a fair comparison."
    CairoMakie.Label(
        fig[3, :],
        disclaimer;
        fontsize = 11,
        color = :gray30,
        word_wrap = true,
        justification = :left,
        halign = :left,
        tellheight = true,
        tellwidth = false,
        padding = (10, 10, 6, 0),
    )

    CairoMakie.save(joinpath(leaderboard_base_path, "boxplot_rmse.png"), fig)
    return nothing
end
