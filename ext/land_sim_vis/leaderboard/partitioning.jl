"""
    PartitionFraction

One row of the partitioning leaderboard: a ratio between sums of diagnostics,
such as `lhf / (lhf + shf)`.

All panels aggregate the numerator and the denominator over the domain or period
of interest and divide only at the end: the mean of a ratio is not the ratio of
the means, and a per-cell `et / precip` diverges as precipitation approaches
zero.
"""
struct PartitionFraction
    """Short name of the ratio, used to build file names."""
    short_name::String

    """Descriptive name shown on the panels."""
    long_name::String

    """The ratio written out, e.g. `"T/ET"`, used to label rows and axes."""
    ratio_label::String

    """Diagnostics summed to form the numerator."""
    numerator::Vector{String}

    """Diagnostics summed to form the denominator."""
    denominator::Vector{String}

    """Components taken from the simulation when forming the *observed* ratio,
    because they are prescribed forcing rather than an independent
    observation."""
    prescribed::Vector{String}

    """Color scale bounds of the bias panel."""
    bias_extrema::Tuple{Float64, Float64}

    """Published global constraint drawn on the seasonal panel when no gridded
    observation exists, as `(value, spread, label)`."""
    reference::Union{Nothing, Tuple{Float64, Float64, String}}
end

"""
The fractions shown in the partitioning leaderboard.

`ef` and `rf` are compared against ERA5. `etp` uses ERA5 evaporation over the
simulation's own precipitation, which is prescribed forcing and so identical on
both sides. `tf` has no gridded counterpart in our artifacts and is shown
against a published global constraint instead.
"""
const PARTITION_FRACTIONS = [
    PartitionFraction(
        "ef",
        "Evaporative Fraction",
        "LHF/(LHF+SHF)",
        ["lhf"],
        ["lhf", "shf"],
        String[],
        (-0.2, 0.2),
        nothing,
    ),
    PartitionFraction(
        "rf",
        "Surface Runoff Fraction",
        "SR/(SR+SSR)",
        ["sr"],
        ["sr", "ssr"],
        String[],
        (-0.3, 0.3),
        nothing,
    ),
    PartitionFraction(
        "etp",
        "Evaporative Index",
        "ET/P",
        ["et"],
        ["precip"],
        ["precip"],
        (-0.3, 0.3),
        nothing,
    ),
    PartitionFraction(
        "tf",
        "Transpiration Fraction",
        "T/ET",
        ["trans"],
        ["et"],
        String[],
        (-0.3, 0.3),
        (0.57, 0.07, "Wei et al. (2017), global annual"),
    ),
]

"""
Cells whose aggregated denominator falls below this fraction of its own global
mean are masked out of the maps: they carry almost no flux, so their ratio is
numerically unstable while meaning nothing physically.
"""
const _MIN_DENOMINATOR_FRACTION = 0.01

"""
    _sum_vars(vars)

Sum a non-empty collection of `OutputVar`s sharing dimensions, returning an
`OutputVar` that keeps the first one's metadata.
"""
function _sum_vars(vars)
    first_var = first(vars)
    length(vars) == 1 && return first_var
    return ClimaAnalysis.remake(first_var; data = sum(var.data for var in vars))
end

"""
    _ratio_var(num, den, fraction)

Build the dimensionless `OutputVar` `num / den`, labelled for `fraction`.
"""
function _ratio_var(num, den, fraction)
    attributes = Dict(
        # `ClimaAnalysis.bias` rejects an empty units string, so dimensionless
        # is spelled out the conventional way rather than left blank.
        "units" => "-",
        "short_name" => fraction.ratio_label,
        "long_name" => fraction.long_name,
    )
    return ClimaAnalysis.remake(num; attributes, data = num.data ./ den.data)
end

"""
    _intersect_mask(base_mask, vars...)

Return a masking function that blanks every cell that `base_mask` removes or
where any of `vars` is missing, so all panels of a row integrate over one common
domain.

The footprint is resolved once and then applied by assignment: resampling the
land-sea mask on every call, as `apply_oceanmask` does, is far too expensive to
repeat for the hundreds of monthly slices a seasonal cycle averages over.
"""
function _intersect_mask(base_mask, vars...)
    blank =
        reduce((a, b) -> a .| b, (isnan.(base_mask(var).data) for var in vars))
    return function (var)
        new_data = copy(var.data)
        new_data[blank] .= NaN
        return ClimaAnalysis.remake(var; data = new_data)
    end
end

"""
    _small_denominator_mask(den, base_mask)

Return a copy of `den` with cells below `_MIN_DENOMINATOR_FRACTION` of the
global mean (and all non-positive cells) set to `NaN`, so `_intersect_mask`
propagates the exclusion to the whole row.
"""
function _small_denominator_mask(den, base_mask)
    masked = base_mask(den)
    global_mean = ClimaAnalysis.weighted_average_lonlat(masked).data[]
    threshold = _MIN_DENOMINATOR_FRACTION * abs(global_mean)
    new_data = copy(masked.data)
    new_data[.!(new_data .> threshold)] .= NaN
    return ClimaAnalysis.remake(masked; data = new_data)
end

"""
    _global_mean_timeseries(var, mask_fn)

Return the lonlat-weighted global mean of `var` at each of its times.
"""
function _global_mean_timeseries(var, mask_fn)
    return [
        ClimaAnalysis.weighted_average_lonlat(
            mask_fn(ClimaAnalysis.slice(var, time = t)),
        ).data[] for t in ClimaAnalysis.times(var)
    ]
end

"""
    _ratio_global_series(diagnostics, fraction, mask_fn)

Return `(dates, num_global, den_global)`: the date of every output time and the
lonlat-weighted global means of `fraction`'s numerator and denominator there,
summed from the `diagnostics` of one side of the comparison.
"""
function _ratio_global_series(diagnostics, fraction, mask_fn)
    num = _sum_vars([diagnostics[c] for c in fraction.numerator])
    den = _sum_vars([diagnostics[c] for c in fraction.denominator])
    return (
        ClimaAnalysis.dates(num),
        _global_mean_timeseries(num, mask_fn),
        _global_mean_timeseries(den, mask_fn),
    )
end

"""
    _monthly_ratio_climatology(dates, num_global, den_global)

Return `(means, spread)`, 12-element vectors indexed by calendar month, of the
global ratio and its standard deviation across years.

Each year contributes one ratio per month, formed from that month's global mean
numerator and denominator; the climatology then averages those ratios over
years. Months sampled in a single year get a spread of zero.
"""
function _monthly_ratio_climatology(dates, num_global, den_global)
    isempty(dates) && return (fill(NaN, 12), fill(NaN, 12))
    months = Dates.month.(dates)
    ratios = num_global ./ den_global
    means, spread = fill(NaN, 12), fill(NaN, 12)
    for m in 1:12
        vals = filter(isfinite, ratios[findall(==(m), months)])
        isempty(vals) && continue
        means[m] = sum(vals) / length(vals)
        spread[m] = length(vals) > 1 ? Statistics.std(vals) : 0.0
    end
    return (means, spread)
end

"""
    _annual_ratios(dates, sim_num, sim_den, obs_num, obs_den)

Return `(years, sim_annual, obs_annual)`: one global ratio per calendar year for
the simulation and for the observations, each summing its numerator and
denominator over the year before dividing, as every other panel does.

Only years where all twelve months are present and finite on both sides
contribute, since a partial year's ratio is aliased by whichever part of the
seasonal cycle it sampled. Each month is weighted by its length.
"""
function _annual_ratios(dates, sim_num, sim_den, obs_num, obs_den)
    years, sim_annual, obs_annual = Int[], Float64[], Float64[]
    series = (sim_num, sim_den, obs_num, obs_den)
    for year in sort(unique(Dates.year.(dates)))
        idxs = findall(d -> Dates.year(d) == year, dates)
        valid = filter(i -> all(isfinite(s[i]) for s in series), idxs)
        length(unique(Dates.month.(dates[valid]))) == 12 || continue
        weights = Dates.daysinmonth.(dates[valid])
        push!(years, year)
        push!(
            sim_annual,
            sum(weights .* sim_num[valid]) / sum(weights .* sim_den[valid]),
        )
        push!(
            obs_annual,
            sum(weights .* obs_num[valid]) / sum(weights .* obs_den[valid]),
        )
    end
    return (years, sim_annual, obs_annual)
end

"""
    _load_partition_sim_var(sim_dir, short_name, spin_up_months)

Load a simulation diagnostic, drop the spin up, and convert it to the units the
observations use.
"""
function _load_partition_sim_var(sim_dir, short_name, spin_up_months)
    var = get(sim_dir, short_name)
    spinup_cutoff = spin_up_months * 31 * 86400.0
    ClimaAnalysis.times(var)[end] >= spinup_cutoff &&
        (var = ClimaAnalysis.window(var, "time", left = spinup_cutoff))
    return preprocess_sim_var(var)
end

"""
    _prepare_partition_row(sim_dir, data_loader, fraction, spin_up_months)

Assemble everything the four panels of one leaderboard row need, or `nothing`
when the simulation does not provide every component of `fraction`.

The observed ratio is built only when every component that is not prescribed
forcing is available from `data_loader`; otherwise the row is model-only.
"""
function _prepare_partition_row(sim_dir, data_loader, fraction, spin_up_months)
    components = union(fraction.numerator, fraction.denominator)
    available = ClimaAnalysis.available_vars(sim_dir)
    missing_components = setdiff(components, available)
    if !isempty(missing_components)
        @info "Skipping $(fraction.short_name): simulation is missing $(join(sort(collect(missing_components)), ", "))"
        return nothing
    end

    # `Any` valued: windowing returns `OutputVar`s backed by views, whose
    # concrete type differs from the freshly loaded ones stored here first.
    sim = Dict{String, Any}(
        c => _load_partition_sim_var(sim_dir, c, spin_up_months) for
        c in components
    )
    start_date = first(values(sim)).attributes["start_date"]

    obs_components = setdiff(components, fraction.prescribed)
    has_obs = issubset(obs_components, available_vars(data_loader))
    obs = Dict{String, Any}()
    if has_obs
        for c in obs_components
            obs_var = get(data_loader, c)
            ClimaAnalysis.set_reference_date!(obs_var, start_date)
            obs[c] = obs_var
        end
        # Window every field to the period all of them cover, then put the
        # observations on the simulation grid.
        all_times = ClimaAnalysis.times.(
            vcat(collect(values(sim)), collect(values(obs))),
        )
        min_time = maximum(first.(all_times))
        max_time = minimum(last.(all_times))
        for d in (sim, obs), c in keys(d)
            d[c] = ClimaAnalysis.window(
                d[c],
                "time",
                left = min_time,
                right = max_time,
            )
        end
        for c in obs_components
            obs[c] = ClimaAnalysis.resampled_as(obs[c], sim[c])
        end
        # Prescribed forcing is identical on both sides by construction.
        for c in fraction.prescribed
            obs[c] = sim[c]
        end
    end

    annual(d, names) =
        ClimaAnalysis.average_time(_sum_vars([d[c] for c in names]))
    sim_num_ann, sim_den_ann =
        annual(sim, fraction.numerator), annual(sim, fraction.denominator)

    base_mask = ClimaAnalysis.apply_oceanmask
    guard_vars =
        Any[sim_num_ann, _small_denominator_mask(sim_den_ann, base_mask)]
    obs_num_ann = obs_den_ann = nothing
    if has_obs
        obs_num_ann, obs_den_ann =
            annual(obs, fraction.numerator), annual(obs, fraction.denominator)
        push!(
            guard_vars,
            obs_num_ann,
            _small_denominator_mask(obs_den_ann, base_mask),
        )
    end
    mask_fn = _intersect_mask(base_mask, guard_vars...)

    sim_ratio = _ratio_var(sim_num_ann, sim_den_ann, fraction)
    obs_ratio =
        has_obs ? _ratio_var(obs_num_ann, obs_den_ann, fraction) : nothing

    lats, sim_num_zonal = _zonal_means(sim_num_ann, mask_fn)
    _, sim_den_zonal = _zonal_means(sim_den_ann, mask_fn)
    # Unlike the mean, the spread has to come from the per-cell ratios: there is
    # no other sense in which a ratio varies along a latitude band.
    # `_small_denominator_mask` has already dropped the cells where one blows up.
    sim_zonal_std = _zonal_std(sim_ratio, mask_fn)
    obs_zonal = obs_zonal_std = nothing
    if has_obs
        _, obs_num_zonal = _zonal_means(obs_num_ann, mask_fn)
        _, obs_den_zonal = _zonal_means(obs_den_ann, mask_fn)
        obs_zonal = obs_num_zonal ./ obs_den_zonal
        obs_zonal_std = _zonal_std(obs_ratio, mask_fn)
    end

    dates, sim_num_global, sim_den_global =
        _ratio_global_series(sim, fraction, mask_fn)
    sim_monthly, sim_spread =
        _monthly_ratio_climatology(dates, sim_num_global, sim_den_global)
    obs_monthly, obs_spread = fill(NaN, 12), fill(NaN, 12)
    years, sim_annual, obs_annual = Int[], Float64[], Float64[]
    if has_obs
        _, obs_num_global, obs_den_global =
            _ratio_global_series(obs, fraction, mask_fn)
        obs_monthly, obs_spread =
            _monthly_ratio_climatology(dates, obs_num_global, obs_den_global)
        years, sim_annual, obs_annual = _annual_ratios(
            dates,
            sim_num_global,
            sim_den_global,
            obs_num_global,
            obs_den_global,
        )
    end

    return (;
        fraction,
        has_obs,
        mask_fn,
        sim_ratio,
        obs_ratio,
        lats,
        sim_zonal = sim_num_zonal ./ sim_den_zonal,
        sim_zonal_std,
        obs_zonal,
        obs_zonal_std,
        sim_monthly,
        sim_spread,
        obs_monthly,
        obs_spread,
        years,
        sim_annual,
        obs_annual,
        year_range = _year_range(first(values(sim))),
    )
end

"""
    _year_range(var)

Return the years spanned by `var` as a display string, e.g. `"2008-2010"`.
"""
function _year_range(var)
    start_date = Dates.DateTime(var.attributes["start_date"])
    times = ClimaAnalysis.times(var)
    y0 = Dates.year(start_date + Dates.Second(round(Int, first(times))))
    y1 = Dates.year(start_date + Dates.Second(round(Int, last(times))))
    return y0 == y1 ? "$(y0)" : "$(y0)–$(y1)"
end

"""
    _fraction_contour_kwargs(row)

Build `more_kwargs` for the SIM panel: fractions live in `[0, 1]`, so the color
scale is fixed there rather than following the data, with extend arrows for the
values (such as an evaporative index above one) that fall outside.
"""
function _fraction_contour_kwargs(row)
    return Dict(
        :plot => Dict(
            :levels => collect(range(0, 1; length = 21)),
            :extendhigh => :auto,
            :extendlow => :auto,
        ),
        :axis => Dict(
            :title => "$(row.fraction.long_name), mean over $(row.year_range)",
        ),
    )
end

"""
    _limit_to_fraction_range!(ax, profiles...)

Set the value axis of a zonal-profile panel so that `[0, 1]` is always in view
and the tails are clipped to the 2nd-98th percentile of the plotted profiles,
band edges included.

Not every fraction is bounded by one: the evaporative fraction diverges wherever
the sensible heat flux offsets the latent heat flux, which happens routinely at
high latitudes and otherwise compresses the whole profile into a sliver of the
axis.
"""
function _limit_to_fraction_range!(ax, profiles...)
    finite = filter(isfinite, vcat(profiles...))
    isempty(finite) && return nothing
    lo = Statistics.quantile(finite, 0.02)
    hi = Statistics.quantile(finite, 0.98)
    CairoMakie.xlims!(ax, min(lo, 0.0), max(hi, 1.0))
    return nothing
end

"""
    _no_obs_panel!(position, fraction)

Fill a panel that has no gridded observations to compare against with a note
saying so.

An empty axis rather than a bare `Label`: a label demands no height, which would
collapse the whole row next to the map panels.
"""
function _no_obs_panel!(position, fraction)
    ax = CairoMakie.Axis(position)
    CairoMakie.hidedecorations!(ax)
    CairoMakie.hidespines!(ax)
    CairoMakie.text!(
        ax,
        0.5,
        0.5;
        text = "No gridded observations\nfor $(fraction.long_name)",
        fontsize = 24,
        align = (:center, :center),
    )
    return nothing
end

"""
    compute_partitioning_leaderboard(leaderboard_base_path,
                                     diagnostics_folder_path;
                                     spin_up_months = 12,
                                     fractions = PARTITION_FRACTIONS,
                                     era5_to_clima_names = ERA5_PARTITION_TO_CLIMA_NAMES)

Plot a leaderboard of the dimensionless fractions that describe how the model
partitions energy and water: the evaporative fraction, the surface runoff
fraction, the evaporative index, and the transpiration fraction.

The layout matches the annual leaderboards produced by
`compute_seasonal_leaderboard`. Rows are fractions and columns are

- `SIM`, the simulated fraction,
- `ANN`, its bias against observations,
- `LAT`, zonal means of the simulated and observed fractions, banded by the
  spread of the per-cell ratios along each latitude band,
- `MON`, their global seasonal cycle with a band showing the spread across
  years,
- `IAV`, the global annual fractions plotted against each other, one point per
  year.

Observations come from ERA5, loaded under the names `era5_to_clima_names` maps.
Rows whose components the simulation does not output are skipped, and rows with
no gridded observations show the simulation alone alongside a published global
constraint.
"""
function compute_partitioning_leaderboard(
    leaderboard_base_path,
    diagnostics_folder_path;
    spin_up_months = 12,
    fractions = PARTITION_FRACTIONS,
    era5_to_clima_names = ERA5_PARTITION_TO_CLIMA_NAMES,
)
    sim_dir = ClimaAnalysis.SimDir(diagnostics_folder_path)
    data_loader = ERA5DataLoader(; era5_to_clima_names)

    rows = filter(
        !isnothing,
        [
            _prepare_partition_row(
                sim_dir,
                data_loader,
                fraction,
                spin_up_months,
            ) for fraction in fractions
        ],
    )
    if isempty(rows)
        @info "No partitioning fractions could be computed from $diagnostics_folder_path"
        return nothing
    end

    groups = ["SIM", "ANN", "LAT", "MON"]
    any(row -> length(row.years) > _MIN_IAV_YEARS, rows) && push!(groups, "IAV")
    fig = CairoMakie.Figure(; size = (600 * length(groups), 400 * length(rows)))

    # The title and the column labels go in before any globe panel: adding a
    # `Label` to a figure that already holds a `GeoAxis` can send GridLayoutBase
    # into infinite recursion. See MakieOrg/GeoMakie.jl#330.
    titlelayout = CairoMakie.GridLayout(
        fig[-1, 1:length(groups)],
        halign = :center,
        tellwidth = false,
    )
    CairoMakie.Label(
        titlelayout[1, 1],
        "Energy and water partitioning\nexcluding spin up ($spin_up_months months)",
        fontsize = 40,
    )
    for (col_idx, group) in enumerate(groups)
        CairoMakie.Label(
            fig[0, col_idx],
            group,
            tellwidth = false,
            fontsize = 30,
        )
    end

    for (row_idx, row) in enumerate(rows)
        CairoMakie.Label(
            fig[row_idx, 0],
            row.fraction.ratio_label,
            tellheight = false,
            fontsize = 30,
        )
        for (col_idx, group) in enumerate(groups)
            if group == "SIM"
                layout = fig[row_idx, col_idx] = CairoMakie.GridLayout()
                ClimaAnalysis.Visualize.contour2D_on_globe!(
                    layout,
                    row.sim_ratio,
                    mask = row.mask_fn;
                    more_kwargs = _fraction_contour_kwargs(row),
                )
            elseif group == "ANN"
                if !row.has_obs
                    _no_obs_panel!(fig[row_idx, col_idx], row.fraction)
                    continue
                end
                layout = fig[row_idx, col_idx] = CairoMakie.GridLayout()
                ClimaAnalysis.Visualize.plot_bias_on_globe!(
                    layout,
                    row.sim_ratio,
                    row.obs_ratio,
                    cmap_extrema = row.fraction.bias_extrema,
                    mask = row.mask_fn,
                )
            elseif group == "IAV"
                if length(row.years) <= _MIN_IAV_YEARS
                    _no_obs_panel!(fig[row_idx, col_idx], row.fraction)
                    continue
                end
                _plot_interannual!(
                    fig[row_idx, col_idx],
                    row.years,
                    row.sim_annual,
                    row.obs_annual;
                    label = row.fraction.ratio_label,
                    show_title = row_idx == 1,
                )
            elseif group == "LAT"
                ax = CairoMakie.Axis(
                    fig[row_idx, col_idx],
                    xlabel = row.fraction.ratio_label,
                    ylabel = "Latitude (degrees)",
                    yticks = -90:30:90,
                    title = row_idx == 1 ? "Zonal mean of the annual mean" : "",
                )
                CairoMakie.ylims!(ax, -90, 90)
                _limit_to_fraction_range!(
                    ax,
                    row.sim_zonal .- row.sim_zonal_std,
                    row.sim_zonal .+ row.sim_zonal_std,
                    row.has_obs ? row.obs_zonal .- row.obs_zonal_std :
                    Float64[],
                    row.has_obs ? row.obs_zonal .+ row.obs_zonal_std :
                    Float64[],
                )
                row.has_obs && _band_zonal_spread!(
                    ax,
                    row.obs_zonal,
                    row.obs_zonal_std,
                    row.lats,
                    :black,
                )
                _band_zonal_spread!(
                    ax,
                    row.sim_zonal,
                    row.sim_zonal_std,
                    row.lats,
                    :firebrick,
                )
                row.has_obs && CairoMakie.lines!(
                    ax,
                    row.obs_zonal,
                    row.lats;
                    color = :black,
                    linewidth = 4,
                    label = "OBS",
                )
                CairoMakie.lines!(
                    ax,
                    row.sim_zonal,
                    row.lats;
                    color = :firebrick,
                    linewidth = 4,
                    label = "SIM",
                )
                # Anchored left: the profiles run up against the right edge,
                # where a legend would sit on top of them.
                row_idx == 1 && CairoMakie.axislegend(
                    ax,
                    position = :lb,
                    framevisible = false,
                )
            else
                ax = CairoMakie.Axis(
                    fig[row_idx, col_idx],
                    xlabel = "Month",
                    ylabel = row.fraction.ratio_label,
                    xticks = (1:12, first.(Dates.monthname.(1:12), 1)),
                    title = row_idx == 1 ? "Global seasonal cycle" : "",
                )
                _plot_partition_reference!(ax, row.fraction.reference)
                if row.has_obs
                    _band_interannual_spread!(
                        ax,
                        row.obs_monthly,
                        row.obs_spread,
                        :black,
                    )
                end
                _band_interannual_spread!(
                    ax,
                    row.sim_monthly,
                    row.sim_spread,
                    :firebrick,
                )
                row.has_obs && CairoMakie.lines!(
                    ax,
                    1:12,
                    row.obs_monthly;
                    color = :black,
                    linewidth = 4,
                    label = "OBS",
                )
                CairoMakie.lines!(
                    ax,
                    1:12,
                    row.sim_monthly;
                    color = :firebrick,
                    linewidth = 4,
                    label = "SIM",
                )
                CairoMakie.axislegend(ax, position = :rt, framevisible = false)
            end
        end
    end

    # The IAV axis is square, so left to a map-sized column it would sit in the
    # middle of a blank margin.
    iav_col = findfirst(==("IAV"), groups)
    isnothing(iav_col) ||
        CairoMakie.colsize!(fig.layout, iav_col, CairoMakie.Fixed(420))

    CairoMakie.save(
        joinpath(leaderboard_base_path, "Partitioning_leaderboard.png"),
        fig,
    )
    return nothing
end

"""
    _plot_partition_reference!(ax, reference)

Shade a published global constraint `(value, spread, label)` across the seasonal
panel. The constraint is an annual global number, so it is drawn flat and
labelled as such rather than implying a month-by-month target.
"""
_plot_partition_reference!(ax, ::Nothing) = nothing

function _plot_partition_reference!(ax, reference)
    value, spread, label = reference
    CairoMakie.band!(
        ax,
        1:12,
        fill(value - spread, 12),
        fill(value + spread, 12);
        color = (:steelblue, 0.15),
        label = label,
    )
    return nothing
end
