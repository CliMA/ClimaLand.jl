# Energy and carbon RMSE boxplots: compare ClimaLand against an ILAMB land-hist
# cohort. "Other-model" RMSE values are inlined from the ILAMB land-hist
# dashboards (https://www.ilamb.org/land-hist/) and from the LE intercomparison
# in the CliMA-Land RMSE tracker. ClimaLand RMSE is computed at runtime from
# the diagnostics directory with ILAMB's default RMSE definition, so that it is
# comparable with the cohort values.

# Convert ET in mm/day to LE in W/m^2 using the latent heat of vaporization
# (2.45e6 J/kg) and seconds per day.
const _ET_TO_LE = 2.45e6 / 86400

# Panel definitions: sim short_name, benchmarks used by the cohort and by
# ClimaLand, and cohort RMSEs.
const _ENERGY_PANELS = (
    (
        title = "H",
        sim_short_name = "shf",
        data_source = "ERA5",
        others_bench = "FLUXCOM",
        sim_bench = "ERA5",
        others = [
            19.1,
            18.0,
            15.3,
            17.9,
            18.9,
            15.8,
            22.9,
            24.0,
            22.6,
            16.8,
            17.8,
            14.8,
        ],
    ),
    (
        title = "LE",
        sim_short_name = "lhf",
        data_source = "ERA5",
        others_bench = "ERA5",
        sim_bench = "ERA5",
        others = [
            0.529,
            0.789,
            0.600,
            0.491,
            0.657,
            0.621,
            0.756,
            0.618,
            0.779,
            0.805,
            0.569,
            0.702,
            0.494,
            0.563,
            0.622,
            0.671,
            0.566,
            0.561,
            0.650,
            0.620,
            0.616,
        ] .* _ET_TO_LE,
    ),
    (
        title = "SWup",
        sim_short_name = "swu",
        data_source = "ERA5",
        others_bench = "CERESed4.1",
        sim_bench = "ERA5",
        others = [
            12.7,
            12.3,
            13.4,
            15.1,
            13.7,
            14.8,
            17.5,
            15.4,
            17.5,
            12.9,
            12.0,
            13.7,
        ],
    ),
    (
        title = "LWup",
        sim_short_name = "lwu",
        data_source = "ERA5",
        others_bench = "CERESed4.1",
        sim_bench = "ERA5",
        others = [12.1, 12.6, 12.3, 13.3, 10.8, 11.3, 11.5, 10.3, 11.0],
    ),
)

const _CARBON_PANELS = (
    (
        title = "GPP",
        sim_short_name = "gpp",
        data_source = "ILAMB",
        others_bench = "FLUXCOM",
        sim_bench = "FLUXCOM",
        others = [
            1.39,
            1.08,
            1.97,
            1.01,
            1.04,
            1.13,
            0.950,
            0.828,
            0.907,
            1.42,
            1.19,
            1.67,
        ],
    ),
    (
        title = "ER",
        sim_short_name = "er",
        data_source = "ILAMB",
        others_bench = "FLUXCOM",
        sim_bench = "FLUXCOM",
        others = [
            1.78,
            1.38,
            2.01,
            4.07,
            3.03,
            3.87,
            4.87,
            2.98,
            4.12,
            2.01,
            2.19,
            2.07,
        ],
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
    _boxplot_rmse(diagnostics_folder_path, short_name, data_source;
                  spin_up_months = 12)

Compute the global RMSE of ClimaLand `short_name` against the benchmark
indicated by `data_source` (`"ILAMB"` or `"ERA5"`), as ILAMB does (see
`_ilamb_cycle_rmse`), so that it is comparable with the cohort values in the
boxplot. The leaderboard figures keep the RMSE of time means. Spinup is removed and both fields are windowed to their
overlap and resampled onto the simulation grid as in
`compute_seasonal_leaderboard`.

Returns `NaN` if the variable is not in the simulation directory or in the
benchmark.
"""
function _boxplot_rmse(
    diagnostics_folder_path,
    short_name,
    data_source;
    spin_up_months = 12,
)
    sim_dir = ClimaAnalysis.SimDir(diagnostics_folder_path)
    data_loader =
        uppercase(data_source) == "ERA5" ? ERA5DataLoader() : ILAMBDataLoader()
    mask_dict = get_mask_dict(data_loader)
    available = ClimaAnalysis.available_vars(
        ClimaAnalysis.SimDir(diagnostics_folder_path),
    )
    (short_name in available_vars(data_loader) && short_name in available) ||
        return NaN

    sim_var = get(sim_dir, short_name)
    obs_var = get(data_loader, short_name)

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

    mask_fn = mask_dict[short_name](sim_var, obs_var)
    return _ilamb_cycle_rmse(sim_var, obs_var, mask_fn)
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
    xlabel =
        "Others vs $(panel.others_bench) (n=$(length(panel.others)))\n" *
        "ClimaLand vs $(panel.sim_bench)"
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

ClimaLand RMSE follows ILAMB's default definition (see `_ilamb_cycle_rmse`).
Energy panels are RMSE in W m⁻² (ClimaLand vs ERA5); the cohort benchmark
differs per panel, as each panel's label shows. Carbon panels are RMSE in
g m⁻² day⁻¹ (ClimaLand and cohort both vs ILAMB FLUXCOM).
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
        rmse_current[p.sim_short_name] = _boxplot_rmse(
            diagnostics_folder_path,
            p.sim_short_name,
            p.data_source,
        )
        rmse_prev[p.sim_short_name] =
            isnothing(prev_diagnostics_folder_path) ? NaN :
            _boxplot_rmse(
                prev_diagnostics_folder_path,
                p.sim_short_name,
                p.data_source,
            )
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
        "strictly fair: reference datasets and forcings differ across " *
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
