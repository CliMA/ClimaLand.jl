"""
    _percentile_contour_kwargs(var; low_q = 0.05, high_q = 0.95, nlevels = 11)

Build a `more_kwargs` dictionary for `ClimaAnalysis.Visualize.contour2D_on_globe!`
that clips the colormap to `[quantile(data, low_q), quantile(data, high_q)]`
and lets `Makie.contourf!` draw extend-arrows for values outside that range.

Returns an empty `Dict` (i.e. fall back to the plotting defaults) when the
field is all-NaN, has fewer than two finite samples, or has a degenerate
percentile range.
"""
function _percentile_contour_kwargs(
    var;
    low_q = 0.05,
    high_q = 0.95,
    nlevels = 11,
)
    vals = filter(isfinite, vec(var.data))
    isempty(vals) && return Dict{Symbol, Any}()
    lo = Statistics.quantile(vals, low_q)
    hi = Statistics.quantile(vals, high_q)
    lo == hi && return Dict{Symbol, Any}()
    levels = collect(range(lo, hi; length = nlevels))
    return Dict(
        :plot =>
            Dict(:levels => levels, :extendhigh => :auto, :extendlow => :auto),
    )
end

"""
    _global_mean_series(sim, obs, mask_fn = ClimaAnalysis.apply_oceanmask)

Return `(dates, sim_global, obs_global)`: the date of every output time and the
lonlat-weighted global mean of `sim` and `obs` there.

Both fields are restricted to the cells where *both* are finite, so the SIM/OBS
gap here matches the global bias in the ANN column even where obs has gaps that
sim does not.
"""
function _global_mean_series(sim, obs, mask_fn = ClimaAnalysis.apply_oceanmask)
    times = ClimaAnalysis.times(sim)
    dates = ClimaAnalysis.dates(sim)
    sim_global = Float64[]
    obs_global = Float64[]
    for t in times
        sim_t = mask_fn(ClimaAnalysis.slice(sim, time = t))
        obs_t = mask_fn(ClimaAnalysis.slice(obs, time = t))
        nan_either = isnan.(sim_t.data) .| isnan.(obs_t.data)
        sim_t.data[nan_either] .= NaN
        obs_t.data[nan_either] .= NaN
        push!(sim_global, ClimaAnalysis.weighted_average_lonlat(sim_t).data[])
        push!(obs_global, ClimaAnalysis.weighted_average_lonlat(obs_t).data[])
    end
    return (dates, sim_global, obs_global)
end

"""
    _monthly_climatology(dates, sim_global, obs_global)

Return `(sim_monthly, obs_monthly, sim_spread, obs_spread)`, each a 12-element
vector indexed by calendar month. The first two are the climatology of the
global means from `_global_mean_series`; the last two their standard deviation
across years, drawn as the interannual band on the MON panel. Months with no
valid sample are `NaN`, months sampled in a single year get a spread of zero.
"""
function _monthly_climatology(dates, sim_global, obs_global)
    isempty(dates) && return ntuple(_ -> fill(NaN, 12), 4)
    months = Dates.month.(dates)
    out_sim, out_obs = fill(NaN, 12), fill(NaN, 12)
    spread_sim, spread_obs = fill(NaN, 12), fill(NaN, 12)
    for m in 1:12
        idxs = findall(==(m), months)
        isempty(idxs) && continue
        for (global_vals, means, spreads) in (
            (sim_global, out_sim, spread_sim),
            (obs_global, out_obs, spread_obs),
        )
            vals = filter(isfinite, global_vals[idxs])
            isempty(vals) && continue
            means[m] = sum(vals) / length(vals)
            spreads[m] = length(vals) > 1 ? Statistics.std(vals) : 0.0
        end
    end
    return (out_sim, out_obs, spread_sim, spread_obs)
end

"""
    _annual_means(dates, sim_global, obs_global)

Return `(years, sim_annual, obs_annual)`, one entry per calendar year of the
global means from `_global_mean_series`.

Only years covering all twelve months in both series contribute: a partial
year's mean is aliased by whichever part of the seasonal cycle it sampled. Each
month is weighted by its length.
"""
function _annual_means(dates, sim_global, obs_global)
    years, sim_annual, obs_annual = Int[], Float64[], Float64[]
    for year in sort(unique(Dates.year.(dates)))
        idxs = findall(d -> Dates.year(d) == year, dates)
        valid = filter(
            i -> isfinite(sim_global[i]) && isfinite(obs_global[i]),
            idxs,
        )
        length(unique(Dates.month.(dates[valid]))) == 12 || continue
        weights = Dates.daysinmonth.(dates[valid])
        push!(years, year)
        push!(sim_annual, sum(weights .* sim_global[valid]) / sum(weights))
        push!(obs_annual, sum(weights .* obs_global[valid]) / sum(weights))
    end
    return (years, sim_annual, obs_annual)
end

"""
    _band_interannual_spread!(ax, means, spread, color)

Shade `means ± spread` over the 12 calendar months on `ax`. Months whose mean is
missing stay blank, and a month sampled in a single year contributes no width.
"""
function _band_interannual_spread!(ax, means, spread, color)
    half_width = [isfinite(s) ? s : 0.0 for s in spread]
    all(iszero, half_width) && return nothing
    CairoMakie.band!(
        ax,
        1:12,
        means .- half_width,
        means .+ half_width;
        color = (color, 0.2),
    )
    return nothing
end

"""
    _band_zonal_spread!(ax, values, spread, lats, color)

Shade `values ± spread` against latitude on `ax`, whose value axis is the
horizontal one. Fully masked latitude bands interrupt the shading rather than
being bridged across.
"""
function _band_zonal_spread!(ax, values, spread, lats, color)
    valid = isfinite.(values) .& isfinite.(spread)
    i = firstindex(valid)
    while i <= lastindex(valid)
        if !valid[i]
            i += 1
            continue
        end
        j = i
        while j < lastindex(valid) && valid[j + 1]
            j += 1
        end
        if j > i
            band = i:j
            CairoMakie.band!(
                ax,
                CairoMakie.Point2f.(values[band] .- spread[band], lats[band]),
                CairoMakie.Point2f.(values[band] .+ spread[band], lats[band]);
                color = (color, 0.2),
            )
        end
        i = j + 1
    end
    return nothing
end

"""
    _zonal_means(var, mask_fn)

Return `(latitudes, values)` for the zonal (longitudinal) mean of the
time-averaged `var` after applying `mask_fn`. Latitudes whose band is entirely
masked or missing come back as `NaN`, which `Makie` skips when drawing the line.

All cells in a latitude band subtend the same area, so no area weighting is
needed.
"""
function _zonal_means(var, mask_fn)
    masked = mask_fn(var)
    zonal = ClimaAnalysis.average_lon(masked)
    return (ClimaAnalysis.latitudes(zonal), vec(zonal.data))
end

"""
    _zonal_std(var, mask_fn)

Return the standard deviation over longitude of the time-averaged `var` within
each latitude band, skipping masked cells the way `_zonal_means` does.
"""
function _zonal_std(var, mask_fn)
    # A band holds the whole population of cells at that latitude rather than a
    # sample of it, so the variance is not Bessel-corrected.
    variance = ClimaAnalysis.variance_lon(mask_fn(var); corrected = false)
    return sqrt.(vec(variance.data))
end

"""
The IAV column is only drawn when some variable has more than this many complete
years; below that a scatter of annual means says nothing about year-to-year
variability.
"""
const _MIN_IAV_YEARS = 3

"""
    _interannual_stats(obs_annual, sim_annual)

Return `(slope, intercept, r2, amplitude_ratio)` for the least squares fit of
the simulated annual means on the observed ones, and `std(sim) / std(obs)`.

Both are reported because they answer different questions: r² is whether the
model varies in the right years, the ratio whether it varies by the right
amount. A model can do the first while damping every anomaly.
"""
function _interannual_stats(obs_annual, sim_annual)
    length(obs_annual) < 3 && return ntuple(_ -> NaN, 4)
    obs_var = Statistics.var(obs_annual)
    sim_var = Statistics.var(sim_annual)
    (obs_var == 0 || sim_var == 0) && return ntuple(_ -> NaN, 4)
    slope = Statistics.cov(obs_annual, sim_annual) / obs_var
    intercept =
        Statistics.mean(sim_annual) - slope * Statistics.mean(obs_annual)
    r2 = Statistics.cor(obs_annual, sim_annual)^2
    return (slope, intercept, r2, sqrt(sim_var / obs_var))
end

"""Ends of the blue → red ramp the IAV points are colored by."""
const _IAV_FIRST_YEAR_COLOR = :royalblue
const _IAV_LAST_YEAR_COLOR = :firebrick

"""
    _label_first_and_last_year!(ax, years, xs, ys)

Write the first and last year beside their own points, in the colors at the ends
of the ramp those points are drawn with. This does the work of a colorbar
without spending a strip of the panel on it, and marks where in the cloud the
run started and ended.
"""
function _label_first_and_last_year!(ax, years, xs, ys)
    x_mid, y_mid = Statistics.mean(xs), Statistics.mean(ys)
    for (idx, color) in (
        (firstindex(years), _IAV_FIRST_YEAR_COLOR),
        (lastindex(years), _IAV_LAST_YEAR_COLOR),
    )
        # Grow the label horizontally back towards the middle; endpoints sit in
        # a corner, where a label growing outwards would be clipped.
        outward = ys[idx] > y_mid
        CairoMakie.text!(
            ax,
            xs[idx],
            ys[idx];
            text = string(years[idx]),
            offset = (0, outward ? 10 : -10),
            align = (
                xs[idx] > x_mid ? :right : :left,
                outward ? :bottom : :top,
            ),
            color,
            fontsize = 14,
        )
    end
    return nothing
end

"""
    _plot_interannual!(position, years, sim_annual, obs_annual;
                       label, show_title)

Draw the IAV panel: one point per year of the area-weighted global annual mean,
simulated against observed, with the 1:1 line, the least squares fit and the
statistics from `_interannual_stats`.

Points run blue to red with the year, so a drift shows up as a march along the
cloud. The axes share one range, putting the 1:1 line at 45°.
"""
function _plot_interannual!(
    position,
    years,
    sim_annual,
    obs_annual;
    label,
    show_title,
)
    ax = CairoMakie.Axis(
        position,
        xlabel = "OBS $label",
        ylabel = "SIM $label",
        aspect = CairoMakie.AxisAspect(1),
        title = show_title ?
                "Annual global means, one point per year\n(dashed 1:1, black least squares fit)" :
                "",
    )
    # Leave room above and below the cloud for the year labels.
    lo, hi = extrema(vcat(sim_annual, obs_annual))
    pad = hi > lo ? 0.14 * (hi - lo) : 1.0
    limits = [lo - pad, hi + pad]
    CairoMakie.limits!(ax, limits..., limits...)
    CairoMakie.lines!(ax, limits, limits; color = :gray, linestyle = :dash)
    slope, intercept, r2, amplitude_ratio =
        _interannual_stats(obs_annual, sim_annual)
    isfinite(slope) && CairoMakie.lines!(
        ax,
        limits,
        intercept .+ slope .* limits;
        color = :black,
        linewidth = 2,
    )
    CairoMakie.scatter!(
        ax,
        obs_annual,
        sim_annual;
        color = years,
        colormap = CairoMakie.cgrad([
            _IAV_FIRST_YEAR_COLOR,
            _IAV_LAST_YEAR_COLOR,
        ]),
        colorrange = (first(years), last(years)),
        markersize = 14,
        strokewidth = 0.5,
        strokecolor = :black,
    )
    _label_first_and_last_year!(ax, years, obs_annual, sim_annual)
    # Annotate the corner the cloud leaves free: points above the 1:1 line fill
    # the upper left, points below it the lower right.
    above = Statistics.mean(sim_annual .- obs_annual) > 0
    CairoMakie.text!(
        ax,
        above ? 0.97 : 0.03,
        above ? 0.03 : 0.97;
        text = Printf.@sprintf(
            "slope %.2f\nr² %.2f\nσ_sim/σ_obs %.2f\n%d years",
            slope,
            r2,
            amplitude_ratio,
            length(years),
        ),
        space = :relative,
        align = above ? (:right, :bottom) : (:left, :top),
        fontsize = 18,
    )
    return nothing
end

"""
    _nee_diverging_contour_kwargs(var, q; nlevels = 21)

Build `more_kwargs` for `contour2D_on_globe!` that uses a diverging
green → white → red colormap centered at zero, so NEE ≈ 0 reads as white,
strong sinks (NEE ≪ 0) as saturated green, and strong sources (NEE ≫ 0)
as saturated red. The half-range is set to the `q`-th percentile of
`|NEE|` so single-cell outliers don't dominate.
"""
function _nee_diverging_contour_kwargs(var, q; nlevels = 21)
    vals = filter(isfinite, vec(var.data))
    isempty(vals) && return Dict{Symbol, Any}()
    bound = Statistics.quantile(abs.(vals), q)
    bound == 0 && return Dict{Symbol, Any}()
    levels = collect(range(-bound, bound; length = nlevels))
    return Dict(
        :plot => Dict(
            :levels => levels,
            :colormap => CairoMakie.cgrad([:darkgreen, :white, :darkred]),
            :extendhigh => :auto,
            :extendlow => :auto,
        ),
    )
end

"""
    _mask_template(var)

Return a lon-lat field on `var`'s grid whose values are all zero.

`_resolved_mask` needs some field on the grid to find out which cells a mask
removes. Zeros are used rather than `var`'s own values so that what comes back
is the mask's footprint alone, with none of `var`'s missing cells in it.
"""
function _mask_template(var)
    lonlat = ClimaAnalysis.slice(var, time = first(ClimaAnalysis.times(var)))
    return ClimaAnalysis.remake(lonlat; data = zeros(size(lonlat.data)))
end

"""
    _resolved_mask(mask_fn, template)

Return a function that masks any field sharing `template`'s grid exactly as
`mask_fn` does, but that resolves *which* cells are masked only once.

`ClimaAnalysis.apply_oceanmask` and the masks from `make_lonlat_mask` resample
onto the target grid on every call, which dominates the runtime of a leaderboard
that masks every monthly slice several times over. The footprint is the same for
any field on the grid. Fields on another grid fall back to `mask_fn` itself.
"""
function _resolved_mask(mask_fn, template)
    blank = isnan.(mask_fn(template).data)
    return function (var)
        size(var.data) == size(blank) || return mask_fn(var)
        new_data = copy(var.data)
        new_data[blank] .= NaN
        return ClimaAnalysis.remake(var; data = new_data)
    end
end

"""
    _prepare_for_bias(base_mask, sim, obs)

Return `(sim, obs, mask_fn)` so that `ClimaAnalysis.bias` / `global_bias` /
`global_rmse` / `plot_bias_on_globe!` integrate over the intersection of finite
cells, matching `_global_mean_series`. `mask_fn` drops every cell that is `NaN`
in either field, so obs with spatial gaps (GOSIF-GPP and residual-ER over
deserts and ice) are averaged over the same area they are normalized by.

Gaps stay `NaN` rather than zero-filled: `resampled_as` erodes the cells next to
one into `NaN` until it is NaN-aware (ClimaAnalysis.jl#198), which loses less
than biasing those edges towards zero.
"""
function _prepare_for_bias(base_mask, sim, obs)
    sim_masked = base_mask(sim)
    obs_masked = base_mask(obs)
    extra_nan = isnan.(sim_masked.data) .| isnan.(obs_masked.data)
    mask_fn = function (var)
        v = base_mask(var)
        any(extra_nan) || return v
        new_data = copy(v.data)
        new_data[extra_nan] .= NaN
        return ClimaAnalysis.remake(v; data = new_data)
    end
    return sim, obs, mask_fn
end

"""
    _get_data_loader(data_source)

Return the observational data loader for `data_source` (case-insensitive), one
of `"ERA5"`, `"FlagshipCarbonMetrics"`, or `"ILAMB"`. Errors on anything else.
"""
function _get_data_loader(data_source)
    source = uppercase(data_source)
    source == "ERA5" && return ERA5DataLoader()
    source == "FLAGSHIPCARBONMETRICS" &&
        return FlagshipCarbonMetricsDataLoader()
    source == "ILAMB" && return ILAMBDataLoader()
    return error(
        "Unknown data_source \"$data_source\"; expected \"ERA5\", \"FlagshipCarbonMetrics\", or \"ILAMB\".",
    )
end

"""
    compute_monthly_leaderboard(leaderboard_base_path,
                                diagnostics_folder_path,
                                data_source)

Plot the biases and a monthly leaderboard of various variables defined over longitude,
latitude, and time. The observational data is determined by `data_source` and can either be
`ILAMB` or `ERA5`.

The parameter `leaderboard_base_path` is the path to save the leaderboards and bias plots,
and `diagnostics_folder_path` is the path to the simulation data.

Loading and preprocessing simulation data is done by loading `OutputVar`s from a
`SimDir` and passing them through `preprocess_sim_var`. Loading and preprocessing
observational data is done by `ERA5DataLoader` or `ILAMBDataLoader`. The masks for
normalizing the global RMSE and bias are determined by `get_mask_dict`. The ranges
of the bias plots are determined by `get_compare_vars_biases_plot_extrema`. See
the functions defined in data_sources.jl.
"""
function compute_monthly_leaderboard(
    leaderboard_base_path,
    diagnostics_folder_path,
    data_source,
)
    sim_dir = ClimaAnalysis.SimDir(diagnostics_folder_path)
    data_loader = _get_data_loader(data_source)
    mask_dict = get_mask_dict(data_loader)

    compare_vars_biases_plot_extrema = get_compare_vars_biases_plot_extrema()
    short_names = intersect(
        ClimaAnalysis.available_vars(sim_dir),
        available_vars(data_loader),
    )
    issubset(short_names, keys(mask_dict)) ||
        error("Not all variables ($short_names) have a mask $(keys(mask_dict))")

    @info "Error against observations"

    # Set up dict for storing simulation and observational data after processing
    # and for storing the month we are interested in
    sim_obs_comparison_dict = Dict()

    # Print dates for debugging
    var = preprocess_sim_var(get(sim_dir, first(short_names)))
    output_dates =
        Dates.DateTime(var.attributes["start_date"]) .+
        Dates.Second.(ClimaAnalysis.times(var))
    @info "Working with dates:"
    @info output_dates

    for short_name in short_names
        @info short_name
        # Simulation data
        sim_var = get(sim_dir, short_name)

        # Observational data
        obs_var = get(data_loader, short_name)

        # Remove first spin_up_months months from simulation
        spin_up_months = 12
        spinup_cutoff = spin_up_months * 30 * 86400.0
        ClimaAnalysis.times(sim_var)[end] >= spinup_cutoff && (
            sim_var =
                ClimaAnalysis.window(sim_var, "time", left = spinup_cutoff)
        )

        # Preprocess sim var to match conventions of data loaders
        sim_var = preprocess_sim_var(sim_var)

        # Make the relative time the same between observational and simulation data
        ClimaAnalysis.set_reference_date!(
            obs_var,
            sim_var.attributes["start_date"],
        )

        # Get the first valid time and last valid time
        obs_times = ClimaAnalysis.times(obs_var)
        sim_times = ClimaAnalysis.times(sim_var)
        min_time = maximum(first.((obs_times, sim_times)))
        max_time = minimum(last.((obs_times, sim_times)))

        # Window OutputVars to restrain the times to those that are the same between
        # both OutputVars
        sim_var = ClimaAnalysis.window(
            sim_var,
            "time",
            left = min_time,
            right = max_time,
        )
        obs_var = ClimaAnalysis.window(
            obs_var,
            "time",
            left = min_time,
            right = max_time,
        )

        obs_var = ClimaAnalysis.shift_longitude(obs_var, -180.0, 180.0)
        obs_var = ClimaAnalysis.resampled_as(obs_var, sim_var)

        sim_obs_comparison_dict[short_name] = (sim_var, obs_var)
    end

    # Plot monthly comparisons
    for short_name in short_names
        sim_var, obs_var = sim_obs_comparison_dict[short_name]

        # Grab the last 12 months if possible for plotting
        times = ClimaAnalysis.times(sim_var)
        min_idx = max(1, length(times) - 11)
        times = times[min_idx:end]
        monthly_dates =
            Dates.DateTime(sim_var.attributes["start_date"]) .+
            Dates.Second.(times)
        num_times = length(times)
        sim_var = ClimaAnalysis.window(
            sim_var,
            "time",
            left = times[begin],
            right = times[end],
        )
        obs_var = ClimaAnalysis.window(
            obs_var,
            "time",
            left = times[begin],
            right = times[end],
        )

        months_and_years = (
            (Dates.monthname(date), Dates.year(date)) for date in monthly_dates
        )

        fig = CairoMakie.Figure(size = (650 * ceil(num_times / 2), 450 * 2))
        mask = _resolved_mask(
            mask_dict[short_name](sim_var, obs_var),
            _mask_template(sim_var),
        )
        times = vcat(
            times,
            Array{Union{Missing, eltype(times)}}(missing, 12 - num_times),
        )
        times = reshape(times, (2, 6))
        for ((indices, t), (month, year)) in zip(pairs(times), months_and_years)
            layout = fig[Tuple(indices)...] = CairoMakie.GridLayout()
            sim_c, obs_c, mask_c = _prepare_for_bias(
                mask,
                ClimaAnalysis.slice(sim_var, time = t),
                ClimaAnalysis.slice(obs_var, time = t),
            )
            ClimaAnalysis.Visualize.plot_bias_on_globe!(
                layout,
                sim_c,
                obs_c,
                cmap_extrema = compare_vars_biases_plot_extrema[short_name],
                mask = mask_c,
            )
            CairoMakie.Label(
                layout[0, 1],
                month * " $year",
                tellwidth = false,
                fontsize = 30,
            )
        end
        CairoMakie.save(
            joinpath(
                leaderboard_base_path,
                "$(data_source)_$(short_name)_bias_plot.png",
            ),
            fig,
        )
    end

    # Plot month (x-axis) and global bias and global RMSE (y-axis)
    fig = CairoMakie.Figure(size = (250 + 450 * length(short_names), 900))
    fig_rmse_bias = fig[1, 1] = CairoMakie.GridLayout()
    for (col, short_name) in enumerate(short_names)
        sim_var, obs_var = sim_obs_comparison_dict[short_name]
        times = ClimaAnalysis.times(sim_var)
        mask = _resolved_mask(
            mask_dict[short_name](sim_var, obs_var),
            _mask_template(sim_var),
        )

        sim_vec = [
            begin
                sim_c, _, mask_c = _prepare_for_bias(
                    mask,
                    ClimaAnalysis.slice(sim_var, time = t),
                    ClimaAnalysis.slice(obs_var, time = t),
                )
                ClimaAnalysis.weighted_average_lonlat(mask_c(sim_c)).data[]
            end for t in times
        ]
        rmse_vec = [
            begin
                sim_c, obs_c, mask_c = _prepare_for_bias(
                    mask,
                    ClimaAnalysis.slice(sim_var, time = t),
                    ClimaAnalysis.slice(obs_var, time = t),
                )
                ClimaAnalysis.global_rmse(sim_c, obs_c, mask = mask_c)
            end for t in times
        ]
        bias_vec = [
            begin
                sim_c, obs_c, mask_c = _prepare_for_bias(
                    mask,
                    ClimaAnalysis.slice(sim_var, time = t),
                    ClimaAnalysis.slice(obs_var, time = t),
                )
                ClimaAnalysis.global_bias(sim_c, obs_c, mask = mask_c)
            end for t in times
        ]

        ax_sim = CairoMakie.Axis(
            fig_rmse_bias[1, col],
            title = "Global sim lonlat averages for $short_name",
            xlabel = "Month",
            ylabel = "Global lonlat averages ($(ClimaAnalysis.units(sim_var)))",
            xticks = (1:12, Dates.monthabbr.(1:12)),
        )
        ax_rmse = CairoMakie.Axis(
            fig_rmse_bias[2, col],
            title = "Global RMSE for $short_name",
            xlabel = "Month",
            ylabel = "Global RMSE ($(ClimaAnalysis.units(sim_var)))",
            xticks = (1:12, Dates.monthabbr.(1:12)),
        )
        ax_bias = CairoMakie.Axis(
            fig_rmse_bias[3, col],
            title = "Global bias for $short_name",
            xlabel = "Month",
            ylabel = "Global bias ($(ClimaAnalysis.units(sim_var)))",
            xticks = (1:12, Dates.monthabbr.(1:12)),
        )

        # Partition months, rmse_vec, and bias_vec into years
        months = Dates.month.(
            Dates.DateTime(sim_var.attributes["start_date"]) .+
            Dates.Second.(times),
        )
        months_split, sim_vec_split, rmse_vec_split, bias_vec_split =
            partition_by_val(12, months, sim_vec, rmse_vec, bias_vec)

        # Plot each year with the earlier year being more transparent than the later years
        axes = (ax_sim, ax_rmse, ax_bias)
        num_years = length(months_split)
        for (curr_year, (months, sim_vec, rmse_vec, bias_vec)) in enumerate(
            zip(months_split, sim_vec_split, rmse_vec_split, bias_vec_split),
        )
            alpha = curr_year / num_years
            data_vecs = [sim_vec, rmse_vec, bias_vec]

            for (ax, data_vec) in zip(axes, data_vecs)
                CairoMakie.lines!(
                    ax,
                    months,
                    data_vec,
                    alpha = alpha,
                    color = :blue,
                )
            end
        end

        # Compute the average over each of the months
        num_months = 12
        num_years = length(sim_vec_split)
        average_per_months = (
            begin
                season_to_avg = compute_group_averages(months, data_vec)
                [get(season_to_avg, month, NaN) for month in 1:12]
            end for data_vec in (sim_vec, rmse_vec, bias_vec)
        )
        for (ax, data_vec) in zip(axes, average_per_months)
            CairoMakie.lines!(ax, 1:12, data_vec, color = :orange)
        end
    end

    # Add a legend for the meaning of the black line
    blue_line = CairoMakie.LineElement(color = :blue)
    orange_line = CairoMakie.LineElement(color = :orange)
    CairoMakie.Legend(
        fig[1, 2],
        [blue_line, orange_line],
        [
            "More transparent - earlier years\nLess transparent - later years",
            "Average over each month",
        ],
    )

    CairoMakie.save(
        joinpath(
            leaderboard_base_path,
            "$(data_source)_global_rmse_and_bias_graphs.png",
        ),
        fig,
    )
end

"""
    compute_seasonal_leaderboard(leaderboard_base_path,
                                 diagnostics_folder_path,
                                 data_source)

Plot the biases and a seasonal leaderboard of various variables defined over longitude,
latitude, and time. The observational data is determined by `data_source` and can either be
`ILAMB` or `ERA5`.

The parameter `leaderboard_base_path` is the path to save the leaderboards and bias plots,
and `diagnostics_folder_path` is the path to the simulation data.

Loading and preprocessing simulation data is done by loading `OutputVar`s from a
`SimDir` and passing them through `preprocess_sim_var`. Loading and preprocessing
observational data is done by `ERA5DataLoader` or `ILAMBDataLoader`. The masks for
normalizing the global RMSE and bias are determined by `get_mask_dict`. The ranges
of the bias plots are determined by `get_compare_vars_biases_plot_extrema`. See
the functions defined in data_sources.jl.
"""
function compute_seasonal_leaderboard(
    leaderboard_base_path,
    diagnostics_folder_path,
    data_source,
)
    # Get everything we need from data_sources.jl
    sim_dir = ClimaAnalysis.SimDir(diagnostics_folder_path)
    data_loader = _get_data_loader(data_source)
    short_names = intersect(
        ClimaAnalysis.available_vars(sim_dir),
        available_vars(data_loader),
    )

    # Need to initialize mask function
    mask_dict = get_mask_dict(data_loader)
    # Store the mask functions after initialization
    mask_fn_dict = Dict()

    compare_vars_biases_plot_extrema = get_compare_vars_biases_plot_extrema()

    # Set up dict for storing simulation and observational data after processing
    # Map short name to Dict which maps season to tuple of OutputVars
    sim_obs_season_comparison_dict = Dict()
    # Map short name to time series of time averages for each season
    sim_obs_time_avg_over_seasons_comparison_dict = Dict()
    # Map short name to the (sim_var, obs_var) full windowed time series, kept
    # for the metadata that survives collapsing along time below.
    sim_obs_full_dict = Dict()
    # Map short name to the global mean at each output time. The MON and IAV
    # columns both reduce it, and the masked slicing behind it is expensive.
    global_series_dict = Dict()
    seasons = ["ANN", "MAM", "JJA", "SON", "DJF"]

    spin_up_months = 12
    for short_name in short_names
        @info short_name
        # Simulation data
        sim_var = get(sim_dir, short_name)

        # Observational data
        obs_var = get(data_loader, short_name)

        # Make masking function
        mask_fn_dict[short_name] = _resolved_mask(
            mask_dict[short_name](sim_var, obs_var),
            _mask_template(sim_var),
        )

        # Remove first spin_up_months from simulation if possible
        spinup_cutoff = spin_up_months * 31 * 86400.0
        ClimaAnalysis.times(sim_var)[end] >= spinup_cutoff && (
            sim_var =
                ClimaAnalysis.window(sim_var, "time", left = spinup_cutoff)
        )

        # Preprocess sim var to match conventions of data loaders
        sim_var = preprocess_sim_var(sim_var)

        # Make the relative time the same between observational and simulation data
        ClimaAnalysis.set_reference_date!(
            obs_var,
            sim_var.attributes["start_date"],
        )

        # Determine which times can be used
        sim_times = ClimaAnalysis.times(sim_var)
        obs_times = ClimaAnalysis.times(obs_var)
        min_time = maximum(first.((sim_times, obs_times)))
        max_time = minimum(last.((sim_times, obs_times)))

        # Window sim_var and obs_var
        sim_var = ClimaAnalysis.window(
            sim_var,
            "time",
            left = min_time,
            right = max_time,
        )
        obs_var = ClimaAnalysis.window(
            obs_var,
            "time",
            left = min_time,
            right = max_time,
        )

        # Resample
        obs_var = ClimaAnalysis.resampled_as(obs_var, sim_var)
        # Reduce along time before collapsing the vars along it below.
        sim_obs_full_dict[short_name] = (sim_var, obs_var)
        global_series_dict[short_name] =
            _global_mean_series(sim_var, obs_var, mask_fn_dict[short_name])
        sim_var_seasons = (sim_var, ClimaAnalysis.split_by_season(sim_var)...)
        obs_var_seasons = (obs_var, ClimaAnalysis.split_by_season(obs_var)...)

        # Take time average
        sim_var_seasons = (
            !isempty(sim_var) ? ClimaAnalysis.average_time(sim_var) : sim_var for sim_var in sim_var_seasons
        )
        obs_var_seasons = (
            !isempty(obs_var) ? ClimaAnalysis.average_time(obs_var) : obs_var for obs_var in obs_var_seasons
        )

        # Save observation and simulation data
        sim_obs_season_comparison_dict[short_name] = Dict(
            season => (sim_var_s, obs_var_s) for
            (season, sim_var_s, obs_var_s) in
            zip(seasons, sim_var_seasons, obs_var_seasons)
        )

        # Compute time averages across seasons
        obs_var_seasons_over_time =
            ClimaAnalysis.split_by_season_across_time(obs_var)
        sim_var_seasons_over_time =
            ClimaAnalysis.split_by_season_across_time(sim_var)

        # Take time average of each season, it is reasonable to assume that
        # there is no missing months between the first and last points in time
        time_averages_sim_var =
            ClimaAnalysis.average_time.(sim_var_seasons_over_time)
        time_averages_obs_var =
            ClimaAnalysis.average_time.(obs_var_seasons_over_time)

        sim_obs_time_avg_over_seasons_comparison_dict[short_name] =
            (time_averages_sim_var, time_averages_obs_var)
    end

    # Make global seasonal averages over all years (after removing spinup)
    # Rows correspond to short names
    # Cols correspond to seasons
    groups = ["ANN", "MAM", "JJA", "SON", "DJF"]
    fig_bias = CairoMakie.Figure(;
        size = (600 * length(groups), 400 * length(short_names)),
    )
    for (row_idx, short_name) in enumerate(short_names)
        CairoMakie.Label(
            fig_bias[row_idx, 0],
            short_name,
            tellheight = false,
            fontsize = 30,
        )
        for (col_idx, group) in enumerate(groups)
            sim_var, obs_var = sim_obs_season_comparison_dict[short_name][group]
            isempty(sim_var) && break
            layout = fig_bias[row_idx, col_idx] = CairoMakie.GridLayout()
            sim_var.attributes["short_name"] = "mean $(ClimaAnalysis.short_name(sim_var))"
            sim_c, obs_c, mask_c =
                _prepare_for_bias(mask_fn_dict[short_name], sim_var, obs_var)
            ClimaAnalysis.Visualize.plot_bias_on_globe!(
                layout,
                sim_c,
                obs_c,
                cmap_extrema = compare_vars_biases_plot_extrema[short_name],
                mask = mask_c,
            )
        end
    end

    # Plot the labels for the short names
    for (col_idx, group) in enumerate(groups)
        CairoMakie.Label(
            fig_bias[0, col_idx],
            group,
            tellwidth = false,
            fontsize = 30,
        )
    end

    # Add a title
    titlelayout = CairoMakie.GridLayout(
        fig_bias[-1, 1:length(groups)],
        halign = :center,
        tellwidth = false,
    )
    CairoMakie.Label(
        titlelayout[1, 1],
        "Annual and seasonal biases over all years excluding spin up ($spin_up_months months)",
        fontsize = 40,
    )
    CairoMakie.save(
        joinpath(
            leaderboard_base_path,
            "$(data_source)_seasonal_avg_over_all_years.png",
        ),
        fig_bias,
    )


    # Add plot of time average for simulation and bias data excluding spin up
    # Rows correspond to short names
    # Cols correspond to "SIM" and "ANN"
    annual_compare_vars_biases_plot_extrema =
        get_compare_vars_biases_plot_extrema(; annual = true)
    annual_means_dict = Dict(
        short_name => _annual_means(global_series_dict[short_name]...) for
        short_name in short_names
    )
    groups = ["SIM", "ANN", "LAT", "MON"]
    any(
        years -> length(years) > _MIN_IAV_YEARS,
        first.(values(annual_means_dict)),
    ) && push!(groups, "IAV")
    fig_sim_ann = CairoMakie.Figure(;
        size = (600 * length(groups), 400 * length(short_names)),
    )

    # Note: Changing the order things are added to this figure
    # may result in an error depending on the GeoMakie version.
    # See https://github.com/MakieOrg/GeoMakie.jl/issues/330

    # Add a title
    titlelayout = CairoMakie.GridLayout(
        fig_sim_ann[-1, 1:length(groups)],
        halign = :center,
        tellwidth = false,
    )
    CairoMakie.Label(
        titlelayout[1, 1],
        "Time average for simulation and bias data\nexcluding spin up ($spin_up_months months)",
        fontsize = 40,
    )

    # Plot the labels for the short names
    for (col_idx, group) in enumerate(groups)
        CairoMakie.Label(
            fig_sim_ann[0, col_idx],
            group,
            tellwidth = false,
            fontsize = 30,
        )
    end

    for (row_idx, short_name) in enumerate(short_names)
        CairoMakie.Label(
            fig_sim_ann[row_idx, 0],
            short_name,
            tellheight = false,
            fontsize = 30,
        )
        for (col_idx, group) in enumerate(groups)
            if group == "SIM"
                sim_var, _ = sim_obs_season_comparison_dict[short_name]["ANN"]
                isempty(sim_var) && break
                layout = fig_sim_ann[row_idx, col_idx] = CairoMakie.GridLayout()
                # NEE's sign is physical, so it gets a diverging colormap
                # centered at zero: green = sink, red = source.
                more_kwargs =
                    short_name == "nee" ?
                    _nee_diverging_contour_kwargs(sim_var, 0.95) :
                    _percentile_contour_kwargs(sim_var)
                # The auto-generated title leaks seconds-since-start_date, and
                # the clipped colorbar hides the actual range of the field.
                sim_var_full, _ = sim_obs_full_dict[short_name]
                sim_dates = ClimaAnalysis.dates(sim_var_full)
                y0, y1 =
                    Dates.year(first(sim_dates)), Dates.year(last(sim_dates))
                year_str = y0 == y1 ? "$(y0)" : "$(y0)–$(y1)"
                # Strip the `, average within …` that ClimaAnalysis appends to
                # `long_name`; the title already says what it averages over.
                long_name =
                    get(sim_var_full.attributes, "long_name", short_name)
                long_name = String(split(long_name, ", average")[1])
                units_str = ClimaAnalysis.units(sim_var)
                masked_data = mask_fn_dict[short_name](sim_var).data
                finite = filter(isfinite, vec(masked_data))
                panel_title = if isempty(finite)
                    "$long_name, mean over $year_str [$units_str]"
                else
                    lo, hi = extrema(finite)
                    Printf.@sprintf(
                        "%s, mean over %s, range %.3g to %.3g [%s]",
                        long_name,
                        year_str,
                        lo,
                        hi,
                        units_str,
                    )
                end
                more_kwargs[:axis] = Dict(:title => panel_title)
                ClimaAnalysis.Visualize.contour2D_on_globe!(
                    layout,
                    sim_var,
                    mask = mask_fn_dict[short_name];
                    more_kwargs,
                )
            elseif group == "LAT"
                sim_var, obs_var =
                    sim_obs_season_comparison_dict[short_name]["ANN"]
                isempty(sim_var) && break
                # Average both fields over the same cells, as the bias plots do,
                # so the two profiles stay comparable where obs has gaps.
                sim_c, obs_c, mask_c = _prepare_for_bias(
                    mask_fn_dict[short_name],
                    sim_var,
                    obs_var,
                )
                sim_lats, sim_zonal = _zonal_means(sim_c, mask_c)
                _, obs_zonal = _zonal_means(obs_c, mask_c)
                sim_zonal_std = _zonal_std(sim_c, mask_c)
                obs_zonal_std = _zonal_std(obs_c, mask_c)
                ax = CairoMakie.Axis(
                    fig_sim_ann[row_idx, col_idx],
                    xlabel = "$short_name ($(ClimaAnalysis.units(sim_var)))",
                    ylabel = "Latitude (degrees)",
                    yticks = -90:30:90,
                    title = row_idx == 1 ? "Zonal mean of the annual mean" : "",
                )
                CairoMakie.ylims!(ax, -90, 90)
                # Bands give the zonal heterogeneity a SIM/OBS gap sits against.
                _band_zonal_spread!(
                    ax,
                    obs_zonal,
                    obs_zonal_std,
                    sim_lats,
                    :black,
                )
                _band_zonal_spread!(
                    ax,
                    sim_zonal,
                    sim_zonal_std,
                    sim_lats,
                    :firebrick,
                )
                # Latitude on the vertical axis, to read alongside the maps.
                CairoMakie.lines!(
                    ax,
                    obs_zonal,
                    sim_lats;
                    color = :black,
                    linewidth = 4,
                    label = "OBS",
                )
                CairoMakie.lines!(
                    ax,
                    sim_zonal,
                    sim_lats;
                    color = :firebrick,
                    linewidth = 4,
                    label = "SIM",
                )
                row_idx == 1 && CairoMakie.axislegend(
                    ax,
                    position = :rb,
                    framevisible = false,
                )
            elseif group == "IAV"
                sim_var_full, _ = sim_obs_full_dict[short_name]
                years, sim_annual, obs_annual = annual_means_dict[short_name]
                length(years) > _MIN_IAV_YEARS || continue
                _plot_interannual!(
                    fig_sim_ann[row_idx, col_idx],
                    years,
                    sim_annual,
                    obs_annual;
                    label = "$short_name ($(ClimaAnalysis.units(sim_var_full)))",
                    show_title = row_idx == 1,
                )
            elseif group == "MON"
                sim_var_full, obs_var_full = sim_obs_full_dict[short_name]
                isempty(sim_var_full) && break
                sim_monthly, obs_monthly, sim_spread, obs_spread =
                    _monthly_climatology(global_series_dict[short_name]...)
                units_str = ClimaAnalysis.units(sim_var_full)
                ax = CairoMakie.Axis(
                    fig_sim_ann[row_idx, col_idx],
                    xlabel = "Month",
                    ylabel = "$short_name ($units_str)",
                    xticks = (
                        1:12,
                        [
                            "J",
                            "F",
                            "M",
                            "A",
                            "M",
                            "J",
                            "J",
                            "A",
                            "S",
                            "O",
                            "N",
                            "D",
                        ],
                    ),
                )
                # Bands show the spread of the global monthly mean across the
                # simulated years, so the SIM/OBS gap can be read against the
                # interannual variability rather than in isolation.
                _band_interannual_spread!(ax, obs_monthly, obs_spread, :black)
                _band_interannual_spread!(
                    ax,
                    sim_monthly,
                    sim_spread,
                    :firebrick,
                )
                CairoMakie.lines!(
                    ax,
                    1:12,
                    obs_monthly;
                    color = :black,
                    linewidth = 4,
                    label = "OBS",
                )
                CairoMakie.lines!(
                    ax,
                    1:12,
                    sim_monthly;
                    color = :firebrick,
                    linewidth = 4,
                    label = "SIM",
                )
                # First row only: on some variables the curves run through the
                # :rt anchor, and repeating the legend clutters the figure.
                row_idx == 1 && CairoMakie.axislegend(
                    ax,
                    position = :rt,
                    framevisible = false,
                )
            else
                sim_var, obs_var =
                    sim_obs_season_comparison_dict[short_name][group]
                isempty(sim_var) && break
                layout = fig_sim_ann[row_idx, col_idx] = CairoMakie.GridLayout()
                sim_c, obs_c, mask_c = _prepare_for_bias(
                    mask_fn_dict[short_name],
                    sim_var,
                    obs_var,
                )
                ClimaAnalysis.Visualize.plot_bias_on_globe!(
                    layout,
                    sim_c,
                    obs_c,
                    cmap_extrema = annual_compare_vars_biases_plot_extrema[short_name],
                    mask = mask_c,
                )
            end
        end
    end

    # The IAV axis is square, so left to a map-sized column it would sit in the
    # middle of a blank margin.
    iav_col = findfirst(==("IAV"), groups)
    isnothing(iav_col) ||
        CairoMakie.colsize!(fig_sim_ann.layout, iav_col, CairoMakie.Fixed(420))

    CairoMakie.save(
        joinpath(
            leaderboard_base_path,
            "$(data_source)_sim_annual_time_avg_over_all_years.png",
        ),
        fig_sim_ann,
    )

    # Make plot with seasons on x-axis and RMSE and bias on the y-axis
    # Rows correspond to short names
    # Cols correspond to SIM and ANN
    # Set up figure to plot on
    fig = CairoMakie.Figure(size = (450 * length(short_names), 900))
    fig_rmse_bias = fig[1, 1] = CairoMakie.GridLayout()


    # Loop over sim_obs_season_over_time_comparison_dict
    for (col, short_name) in enumerate(short_names)
        sim_vars, obs_vars =
            sim_obs_time_avg_over_seasons_comparison_dict[short_name]
        mask = mask_fn_dict[short_name]

        # Get season and compute global bias and global rmse
        seasons = [sim_var.attributes["season"] for sim_var in sim_vars]
        sim_vec = [
            begin
                sim_c, _, mask_c = _prepare_for_bias(mask, sim_var, obs_var)
                ClimaAnalysis.weighted_average_lonlat(mask_c(sim_c)).data[]
            end for (sim_var, obs_var) in zip(sim_vars, obs_vars)
        ]
        rmse_vec = [
            begin
                sim_c, obs_c, mask_c =
                    _prepare_for_bias(mask, sim_var, obs_var)
                ClimaAnalysis.global_rmse(sim_c, obs_c, mask = mask_c)
            end for (sim_var, obs_var) in zip(sim_vars, obs_vars)
        ]
        bias_vec = [
            begin
                sim_c, obs_c, mask_c =
                    _prepare_for_bias(mask, sim_var, obs_var)
                ClimaAnalysis.global_bias(sim_c, obs_c, mask = mask_c)
            end for (sim_var, obs_var) in zip(sim_vars, obs_vars)
        ]

        # Partition by seasons
        # Map each season to a number for plotting
        season_to_num = Dict("MAM" => 1, "JJA" => 2, "SON" => 3, "DJF" => 4)
        seasons = [season_to_num[season] for season in seasons]
        seasons_split, sim_vec_split, bias_vec_split, rmse_vec_split =
            partition_by_val(4, seasons, sim_vec, bias_vec, rmse_vec)

        # Set up three axes
        ax_sim = CairoMakie.Axis(
            fig_rmse_bias[1, col],
            title = "Global sim lonlat averages for $short_name",
            xlabel = "Season",
            ylabel = "Global lonlat averages ($(ClimaAnalysis.units(first(sim_vars))))",
            xticks = (1:4, ["MAM", "JJA", "SON", "DJF"]),
        )
        ax_rmse = CairoMakie.Axis(
            fig_rmse_bias[2, col],
            title = "Global RMSE for $short_name",
            xlabel = "Season",
            ylabel = "Global RMSE ($(ClimaAnalysis.units(first(sim_vars))))",
            xticks = (1:4, ["MAM", "JJA", "SON", "DJF"]),
        )
        ax_bias = CairoMakie.Axis(
            fig_rmse_bias[3, col],
            title = "Global bias for $short_name",
            xlabel = "Season",
            ylabel = "Global bias ($(ClimaAnalysis.units(first(sim_vars))))",
            xticks = (1:4, ["MAM", "JJA", "SON", "DJF"]),
        )

        # Plot on axes
        axes = (ax_sim, ax_rmse, ax_bias)
        num_years = length(seasons_split)
        for (curr_year, (seasons, sim_vec, rmse_vec, bias_vec)) in enumerate(
            zip(seasons_split, sim_vec_split, rmse_vec_split, bias_vec_split),
        )
            alpha = curr_year / num_years
            data_vecs = [sim_vec, rmse_vec, bias_vec]

            for (ax, data_vec) in zip(axes, data_vecs)
                CairoMakie.lines!(
                    ax,
                    seasons,
                    data_vec,
                    alpha = alpha,
                    color = :blue,
                )
            end
        end

        # Compute the average over each of the seasons
        num_years = length(sim_vec_split)
        average_per_seasons = (
            begin
                season_to_avg = compute_group_averages(seasons, data_vec)
                [get(season_to_avg, season, NaN) for season in 1:4]
            end for data_vec in (sim_vec, rmse_vec, bias_vec)
        )

        for (ax, data_vec) in zip(axes, average_per_seasons)
            CairoMakie.lines!(ax, 1:4, data_vec, color = :orange)
        end
    end

    # Add a legend for the meaning of the black line
    blue_line = CairoMakie.LineElement(color = :blue)
    orange_line = CairoMakie.LineElement(color = :orange)
    CairoMakie.Legend(
        fig[1, 2],
        [blue_line, orange_line],
        [
            "More transparent - earlier years\nLess transparent - later years",
            "Average over each season",
        ],
    )

    CairoMakie.save(
        joinpath(
            leaderboard_base_path,
            "$(data_source)_seasonal_global_rmse_and_bias_graphs.png",
        ),
        fig,
    )
end

"""
    partition_by_val(val, x_vec, y_vecs...)

Partition `x_vec` and `y_vecs` into subarrays, so that the last value of each
subarray is `val` if possible.

Examples
=========

```jldoctest
julia> partition_by_val(4, [3,4,1,2,3,4,1,2], collect(1:8))
2-element Vector{Vector{Vector{Int64}}}:
 [[3, 4], [1, 2, 3, 4], [1, 2]]
 [[1, 2], [3, 4, 5, 6], [7, 8]]
```
"""
function partition_by_val(val, x_vec, y_vecs...)
    for vec in y_vecs
        length(x_vec) != length(vec) && error(
            "Length of $x_vec ($(length(x_vec))) and length of $vec ($(length(vec))) are not the same",
        )
    end
    start_and_end_indices = vcat(0, findall(==(val), x_vec))
    if !(length(x_vec) in start_and_end_indices)
        push!(start_and_end_indices, length(x_vec))
    end
    vecs = (x_vec, y_vecs...)
    ret_vecs = [Vector{eltype(vec)}[] for vec in vecs]
    for idx in eachindex(start_and_end_indices[1:(end - 1)])
        for (vec, ret_vec) in zip(vecs, ret_vecs)
            start_idx = start_and_end_indices[idx] + 1
            end_idx = start_and_end_indices[idx + 1]
            push!(ret_vec, vec[start_idx:end_idx])
        end
    end
    return ret_vecs
end

"""
    compute_group_averages(groups, vals)

Return a dictionary mapping values in `groups` to average of the corresponding
values in `vals`.

Examples
=========

```jldoctest
julia> compute_group_averages([1,2,3,4,1,2], [1,2,100,200,11,-2])
Dict{Int64, Float64} with 4 entries:
  4 => 200.0
  2 => 0.0
  3 => 100.0
  1 => 6.0
```
"""
function compute_group_averages(groups, vals)
    length(groups) != length(vals) &&
        error("Length of $groups and $vals are not the same")
    group_to_vals = Dict{eltype(groups), Vector{eltype(vals)}}()
    for (group, val) in zip(groups, vals)
        if group ∉ keys(group_to_vals)
            group_to_vals[group] = [val]
        else
            push!(group_to_vals[group], val)
        end
    end
    return Dict(group => mean(vals) for (group, vals) in group_to_vals)
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 3
        error(
            "Usage: julia leaderboard.jl <Filepath to save leaderboard and bias plots> <Filepath to simulation data> <\"ERA5\" or \"ILAMB\">",
        )
    end
    leaderboard_base_path = ARGS[begin]
    diagnostics_folder_path = ARGS[begin + 1]
    data_source = ARGS[begin + 2]

    # Error handling
    isdir(leaderboard_base_path) ||
        error("$leaderboard_base_path is not a directory")

    compute_monthly_leaderboard(
        leaderboard_base_path,
        diagnostics_folder_path,
        data_source,
    )
    compute_seasonal_leaderboard(
        leaderboard_base_path,
        diagnostics_folder_path,
        data_source,
    )
end
