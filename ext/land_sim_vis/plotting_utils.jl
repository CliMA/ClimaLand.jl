"""
    make_heatmaps(
        savedir,
        diagdir,
        short_names,
        date;
        plot_stem_name = "figures",
        levels = nothing,
        plot! = viz.heatmap2D_on_globe!,
        mask = viz.oceanmask(),
        plot_kwargs = Dict(
            :mask => ClimaAnalysis.Utils.kwargs(color = :white),
            :plot => ClimaAnalysis.Utils.kwargs(rasterize = true),
        ),
    )

Generates png files called in the provided savedir;
these images are global maps for the provided short_names variables
contained in the provided savedir folder (see ClimaAnalysis documentation),
at the date provided. Plots are saved under the name
short_name_date_plot_stem_name.png.
Variables with multiple levels are plotted at each level
in levels, with the level appended to the plot name;
if levels is nothing, the surface is plotted by default.

The plotting defaults are set for global plots with an ocean mask.
"""
function LandSimVis.make_heatmaps(
    savedir,
    diagdir,
    short_names,
    date;
    plot_stem_name = "figures",
    levels = nothing,
    plot! = viz.heatmap2D_on_globe!,
    mask = viz.oceanmask(),
    plot_kwargs = Dict(
        :mask => ClimaAnalysis.Utils.kwargs(color = :white),
        :plot => ClimaAnalysis.Utils.kwargs(rasterize = true),
    ),
)
    simdir = ClimaAnalysis.SimDir(diagdir)
    for short_name in short_names
        var = get(simdir; short_name)
        avail_dates = ClimaAnalysis.dates(var)
        idx = argmin(abs.(date .- avail_dates))
        if ClimaAnalysis.has_altitude(var)
            plot_zvals =
                levels isa Nothing ? [ClimaAnalysis.altitudes(var)[end]] :
                ClimaAnalysis.altitudes(var)[levels] # sfc or levels chosen
            for zval in plot_zvals
                kwarg_z = Dict(:z => zval)
                fig = CairoMakie.Figure(size = (600, 400))
                if mask isa Nothing
                    plot!(
                        fig,
                        ClimaAnalysis.slice(
                            var,
                            time = ClimaAnalysis.times(var)[idx];
                            kwarg_z...,
                        ),
                        more_kwargs = plot_kwargs,
                    )
                else
                    plot!(
                        fig,
                        ClimaAnalysis.slice(
                            var,
                            time = ClimaAnalysis.times(var)[idx];
                            kwarg_z...,
                        ),
                        mask = mask,
                        more_kwargs = plot_kwargs,
                    )
                end
                CairoMakie.save(
                    joinpath(
                        savedir,
                        "$(short_name)_$(zval)_$(avail_dates[idx])_$(plot_stem_name).png",
                    ),
                    fig,
                )
            end
        else
            kwarg_z = Dict()
            fig = CairoMakie.Figure(size = (600, 400))
            if mask isa Nothing
                plot!(
                    fig,
                    ClimaAnalysis.slice(
                        var,
                        time = ClimaAnalysis.times(var)[idx];
                        kwarg_z...,
                    ),
                    more_kwargs = plot_kwargs,
                )
            else
                plot!(
                    fig,
                    ClimaAnalysis.slice(
                        var,
                        time = ClimaAnalysis.times(var)[idx];
                        kwarg_z...,
                    ),
                    more_kwargs = plot_kwargs,
                    mask = mask,
                )
            end

            CairoMakie.save(
                joinpath(
                    savedir,
                    "$(short_name)_$(avail_dates[idx])_$(plot_stem_name).png",
                ),
                fig,
            )
        end
    end
    return nothing
end

"""
    check_conservation(savedir, diagdir; plot_stem_name = "conservation")

Generates two png files assessing water and energy conservation in
savedir/water_plot_stem_name.png
savedir/energy_plot_stem_name.png

The resulting png files contain a time series of the global mean
(area-weighted) energy or water volume error, in units of
`J/m^2` and `m`. Only continents are included in the global average.
"""
function LandSimVis.check_conservation(
    savedir,
    diagdir;
    plot_stem_name = "conservation",
)
    simdir = ClimaAnalysis.SimDir(diagdir)
    ## Energy
    energy_per_area = get(simdir; short_name = "epa")
    energy_per_area_change = get(simdir; short_name = "epac")
    N = length(ClimaAnalysis.times(energy_per_area)) - 1
    times = ClimaAnalysis.times(energy_per_area)
    energy_0 = ClimaAnalysis.apply_oceanmask(
        ClimaAnalysis.slice(energy_per_area; time = times[1]),
    )
    energy_end = ClimaAnalysis.apply_oceanmask(
        ClimaAnalysis.slice(energy_per_area; time = times[end]),
    )
    mean_energy =
        (
            ClimaAnalysis.weighted_average_lonlat(energy_0).data[1] +
            ClimaAnalysis.weighted_average_lonlat(energy_end).data[1]
        ) / 2

    water_volume_per_area = get(simdir; short_name = "wvpa")
    water_volume_per_area_change = get(simdir; short_name = "wvpac")
    water_volume_0 = ClimaAnalysis.apply_oceanmask(
        ClimaAnalysis.slice(water_volume_per_area; time = times[1]),
    )
    water_volume_end = ClimaAnalysis.apply_oceanmask(
        ClimaAnalysis.slice(water_volume_per_area; time = times[end]),
    )
    mean_water_volume =
        (
            ClimaAnalysis.weighted_average_lonlat(water_volume_0).data[1] +
            ClimaAnalysis.weighted_average_lonlat(water_volume_end).data[1]
        ) / 2
    energy_error = zeros(N)
    water_volume_error = zeros(N)
    for i in 1:1:N
        # error = nanmean[(X(t) - X(0) - Expected Change in X)]
        t = times[i]
        #tp1 = times[i + 1]
        energy_error[i] =
            ClimaAnalysis.weighted_average_lonlat(
                ClimaAnalysis.apply_oceanmask(
                    ClimaAnalysis.slice(energy_per_area, time = t),
                ) - energy_0 - ClimaAnalysis.apply_oceanmask(
                    ClimaAnalysis.slice(energy_per_area_change, time = t),
                ),
            ).data[1] ./ mean_energy

        water_volume_error[i] =
            ClimaAnalysis.weighted_average_lonlat(
                ClimaAnalysis.apply_oceanmask(
                    ClimaAnalysis.slice(water_volume_per_area, time = t),
                ) - water_volume_0 - ClimaAnalysis.apply_oceanmask(
                    ClimaAnalysis.slice(
                        water_volume_per_area_change,
                        time = t,
                    ),
                ),
            ).data[1] ./ mean_water_volume

    end
    titles =
        ["Global mean energy per area", "Global mean water volume per area"]
    quantity_names = ["energy", "water"]
    errors = [energy_error, water_volume_error]
    typical_value =
        [@sprintf("%1.2le", mean_energy), @sprintf("%1.2le", mean_water_volume)]
    units = ["J/m²", "m³/m²"]
    for i in 1:2
        fig_cycle = CairoMakie.Figure(size = (600, 400))
        ax = Axis(
            fig_cycle[1, 1],
            xlabel = "Years",
            ylabel = "Global Mean Fractional Conservation Error",
            title = "$(titles[i]), typical value = $(typical_value[i]) $(units[i])",
        )
        CairoMakie.lines!(ax, times[1:N] ./ 24 ./ 3600 ./ 365, errors[i])
        CairoMakie.save(
            joinpath(savedir, "$(quantity_names[i])_$(plot_stem_name).png"),
            fig_cycle,
        )
    end
    return nothing
end

"""
    make_ocean_masked_annual_timeseries(
        savedir,
        diagdir,
        short_names;
	plot_stem_name = "annual_timeseries"
    )

Generates multiple png files called short_name_plot_stem_name.png in the provided savedir,
one for each short_name in the Vector provided `short_names`.
These images contain the timeseries for the global mean of the provided short_names variables
contained in the provided savedir folder (see ClimaAnalysis documentation).
"""
function make_ocean_masked_annual_timeseries(
    savedir,
    diagdir,
    short_names;
    plot_stem_name = "annual_timeseries",
)
    simdir = ClimaAnalysis.SimDir(diagdir)
    for short_name in short_names
        var = get(simdir; short_name)
        kwarg_z = ClimaAnalysis.has_altitude(var) ? Dict(:z => 1) : Dict() # if has altitude, take first layer
        var_sliced = ClimaAnalysis.slice(var; kwarg_z...)
        # var_global_average below is a vector of vector, one for each year of simulation, containing monthly global average of var.
        # i represent a year, from 1 to last year
        # for more details on the ClimaAnalysis functions, see ClimaAnalysis docs.
        var_global_average = [
            ClimaAnalysis.weighted_average_lonlat(
                ClimaAnalysis.apply_oceanmask(
                    ClimaAnalysis.window(
                        var_sliced,
                        "time",
                        left = (i - 1) * 366 * 86400 + 30 * 86400, # 1 year left of year i, in seconds.
                        right = i * 366 * 86400, # 1 year right of year i, in seconds
                    ),
                ),
            ).data for i in range(
                1,
                round(last(ClimaAnalysis.times(var_sliced) / (365 * 86400))), # n years
            )
        ]
        fig_seasonal_cycle = CairoMakie.Figure(size = (600, 400))
        ax = Axis(
            fig_seasonal_cycle[1, 1],
            xlabel = "Month of the year",
            ylabel = "$short_name [$(ClimaAnalysis.units(var))]",
            xticks = (1:1:12, Dates.monthabbr.(1:12)),
        )
        [
            CairoMakie.lines!(
                ax,
                var_global_average[i],
                color = RGBf(0.5, 0.5, 0.5),
                linewidth = 1,
                linestyle = (i == 1 ? :dash : :solid), # dashed line for the 1st year
            ) for i in 1:length(var_global_average)
        ]
        # The three next lines here are computing the average for each month of var, January to December. It accounts for cases if the last simulated year is incomplete. In that case, the last vector of var_global_average would be shorter than 12, so in order to compute the average, it needs to be padded with NaNs until it reaches the length of 12.
        # We also want the average for each month to be computed even if there is NaNs in the last vector, for example January may have 3 data points, but December only two, but we still want the average of all months.
        max_len = maximum(length.(var_global_average)) # 12 months
        padded_var = [
            vcat(v, fill(NaN, max_len - length(v))) for v in var_global_average
        ] # fill with NaN up to 12, if a vector is shorter
        monthly_averages_var = [
            mean(filter(!isnan, collect([v[i] for v in padded_var]))) for
            i in 1:max_len
        ] # compute average of each month, ignoring potential NaNs
        CairoMakie.lines!(
            ax,
            monthly_averages_var,
            color = :black,
            linewidth = 3,
        )
        CairoMakie.xlims!(ax, 1, 12)
        CairoMakie.save(
            joinpath(savedir, "$(short_name)_$(plot_stem_name).png"),
            fig_seasonal_cycle,
        )
    end
    return nothing
end

"""
   time_to_date(t::AbstractFloat, start_date)

Converts a time since the start_date (measured in seconds)
to a date.
"""
function time_to_date(t::AbstractFloat, start_date)
    return start_date + Dates.Millisecond(round(1_000 * t))
end

"""
   time_to_date(t::ITime, start_date)

Converts an ITime to a date using the epoch
of the Itime, the counter, and the period (unit)
of the counter.

Although the epoch can be different from the start_date,
we usually think of the simulation time as relative to the start_date,
and so we warn here if that is not the case
"""
function time_to_date(t::ITime, start_date)
    start_date != t.epoch &&
        @warn("$(start_date) is different from the simulation time epoch.")
    return isnothing(t.epoch) ? start_date + t.counter * t.period : date(t)
end

"""
    _reduction_starts(dates, start_date)

Move the timestamps of a reduction (a monthly average, say) from the end of each
averaging period to its start, so that output kept in memory is dated the way
output written to disk is.

The two writers disagree on when a reduction happened. A `DictWriter` stores it
under the time the period closed, while `ClimaDiagnostics` writes the start of
the period to NetCDF: `init_time` for the first output, the end of the previous
period after that. Reading the two paths through `site_timeseries` therefore
needs this on the `DictWriter` side, and only for reductions — an instantaneous
diagnostic carries no period and is already stamped at its own time.
"""
function _reduction_starts(dates, start_date)
    isempty(dates) && return dates
    return [start_date; dates[1:(end - 1)]]
end

"""
    site_source(diagnostics)

Return what the diagnostics of a single-site simulation are read from: a
`ClimaAnalysis.SimDir` over the output directory when they were written to disk,
or the `DictWriter` holding them when they were kept in memory. As everywhere
else here, the whole list is taken to share one writer.

Build this once and hand it to `site_timeseries` when reading several variables.
A `SimDir` walks the output tree and then caches every `OutputVar` it reads, so
rebuilding it per variable does that walk again and throws the cache away.
"""
function site_source(diagnostics)
    writer = first(diagnostics).output_writer
    writer isa ClimaDiagnostics.Writers.NetCDFWriter || return writer
    return ClimaAnalysis.SimDir(writer.output_dir)
end

"""
    site_timeseries(diagnostics, short_name, start_date;
                    layer = nothing, source = site_source(diagnostics))

Return `(dates, values, units)` for the single-site diagnostic `short_name`,
whether the simulation stored its output in memory with a `DictWriter` or on
disk with a `NetCDFWriter`.

Dates label the start of each reduction period in both cases, so plots do not
have to know which writer produced them.

For variables resolved in depth, `layer` selects a level counting from `1` at
the bottom; the default is the top level.
"""
function site_timeseries(
    diagnostics,
    short_name,
    start_date;
    layer = nothing,
    source = site_source(diagnostics),
)
    return _site_timeseries(source, diagnostics, short_name, start_date, layer)
end

# On disk the SimDir holds everything: it knows the units, and the dates it
# reports are already the start of each reduction period.
function _site_timeseries(
    sim_dir::ClimaAnalysis.SimDir,
    _,
    short_name,
    _,
    layer,
)
    var = get(sim_dir, short_name)
    if ClimaAnalysis.has_altitude(var)
        altitudes = ClimaAnalysis.altitudes(var)
        layer_id = isnothing(layer) ? length(altitudes) : layer
        var = ClimaAnalysis.slice(var, z = altitudes[layer_id])
    end
    return (ClimaAnalysis.dates(var), vec(var.data), ClimaAnalysis.units(var))
end

# In memory there is no file to describe the output, so the diagnostic itself
# supplies the name it is stored under, its units, and whether its timestamps
# need `_reduction_starts`.
function _site_timeseries(
    writer::ClimaDiagnostics.Writers.DictWriter,
    diagnostics,
    short_name,
    start_date,
    layer,
)
    matches = filter(d -> d.variable.short_name == short_name, diagnostics)
    isempty(matches) &&
        error("$short_name is not available in the saved diagnostics.")
    diagnostic = first(matches)
    times, values = ClimaLand.Diagnostics.diagnostic_as_vectors(
        writer,
        diagnostic.output_short_name;
        layer,
    )
    dates = time_to_date.(times, start_date)
    isnothing(diagnostic.reduction_time_func) ||
        (dates = _reduction_starts(dates, start_date))
    return (dates, values, diagnostic.variable.units)
end

"""
    make_diurnal_timeseries(
    savedir,
    diagnostics,
    start_date;
    plot_stem_name = "diurnal_timeseries",
    comparison_data=nothing,
    spinup_date=start_date,
)

Creates multiple png files, each showing the
average diurnal cycle for the diagnostics; the files
is saved under short_name_plot_stem_name.png in the directory `savedir`.

The `start_date` is used to convert from timestamps of seconds since
the start date (in diagnostics) to dates; only values observed or
simulated after the `spinup_date` are included.

To include a comparison to data, a NamedTuple `comparison_data`
may optionally be passed. This should include the timeseries of the
data labeled with the same variable name as the diagnostics use. For example:
comparison_data = (; UTC_datetime, gpp = [....], shf = [....]) will
result in the timeseries of gpp vs UTC_datetime and shf vs UTC_datetime
being plotted, provided that those diagnostics were recorded during the simulation.
"""
function LandSimVis.make_diurnal_timeseries(
    savedir,
    diagnostics,
    start_date;
    plot_stem_name = "diurnal_timeseries",
    comparison_data = nothing,
    spinup_date = start_date,
)
    source = site_source(diagnostics)
    for d in diagnostics
        sn = d.variable.short_name
        model_dates, model_output, unit =
            site_timeseries(diagnostics, sn, start_date; source)
        spinup_idx = findfirst(spinup_date .<= model_dates)
        hour_of_day, model_diurnal_cycle = compute_diurnal_cycle(
            model_dates[spinup_idx:end],
            model_output[spinup_idx:end],
        )
        fig = CairoMakie.Figure(size = (800, 400))
        ax = CairoMakie.Axis(
            fig[1, 1],
            xlabel = "Hour of day (UTC)",
            ylabel = "$sn [$(unit)]",
        )
        CairoMakie.lines!(
            ax,
            hour_of_day,
            model_diurnal_cycle,
            label = "Model",
            color = "blue",
        )
        if ~(comparison_data isa Nothing) &&
           (Symbol(sn) ∈ propertynames(comparison_data))
            data = getproperty(comparison_data, Symbol(sn))
            data_dates = getproperty(comparison_data, :UTC_datetime)
            spinup_idx = findfirst(spinup_date .<= data_dates)
            hour_of_day, data_diurnal_cycle = compute_diurnal_cycle(
                data_dates[spinup_idx:end],
                data[spinup_idx:end],
            )
            RMSD = @sprintf "%.2e" sqrt(
                mean(abs2, model_diurnal_cycle .- data_diurnal_cycle),
            )

            R² = cor(model_diurnal_cycle, data_diurnal_cycle)^2
            CairoMakie.lines!(
                ax,
                hour_of_day,
                data_diurnal_cycle,
                label = "Data",
                color = "yellow",
            )
            ax.title =
                "$(sn): RMSD = " *
                RMSD *
                ", R² = $(round(R²[1][1], digits = 2))"
        end
        axislegend(ax, position = :lt)
        CairoMakie.save(joinpath(savedir, "$(sn)_$(plot_stem_name).png"), fig)
    end
    return nothing
end

"""
    compute_diurnal_cycle(dates, data)

Computes the average diurnal cycle by binning the data dates into hour
intervals and then computing the mean value of the data per interval.
Returns the hours of the day and the mean value of the data each hour.
"""
function compute_diurnal_cycle(dates, data)
    hour_of_day = Hour.(dates)
    mean_by_hour = [mean(data[hour_of_day .== Hour(i)]) for i in 0:23]
    return [Hour(i) for i in 0:23], mean_by_hour
end

"""
    make_timeseries(
    savedir,
    diagnostics,
    start_date;
    layer = nothing, 
    plot_stem_name = "timeseries",
    comparison_data=nothing,
    spinup_date=start_date,
)

Creates multiple png files, each showing the
timeseries for the diagnostics; the files
is saved under short_name_plot_stem_name.png in the directory `savedir`.

The `start_date` is used to convert from timestamps of seconds since
the start date (in diagnostics) to dates; only values observed or
simulated after the `spinup_date` are included.

To include a comparison to data, a NamedTuple `comparison_data`
may optionally be passed. This should include the timeseries of the
data labeled with the same variable name as the diagnostics use. For example:
comparison_data = (; UTC_datetime, gpp = [....], shf = [....]) will
result in the timeseries of gpp vs UTC_datetime and shf vs UTC_datetime
being plotted, provided that those diagnostics were recorded during the simulation.
"""
function LandSimVis.make_timeseries(
    savedir,
    diagnostics,
    start_date;
    layer = nothing,
    plot_stem_name = "timeseries",
    comparison_data = nothing,
    spinup_date = start_date,
)
    source = site_source(diagnostics)
    for d in diagnostics
        sn = d.variable.short_name
        model_dates, model_output, unit =
            site_timeseries(diagnostics, sn, start_date; layer, source)
        spinup_idx = findfirst(spinup_date .<= model_dates)
        fig = CairoMakie.Figure(size = (800, 400))
        ax = CairoMakie.Axis(
            fig[1, 1],
            xlabel = "Date (UTC)",
            ylabel = "$sn [$(unit)]",
        )
        CairoMakie.lines!(
            ax,
            model_dates[spinup_idx:end],
            model_output[spinup_idx:end],
            label = "Model",
            color = "blue",
        )
        xlims = extrema(model_dates[spinup_idx:end])
        xlims!(ax, xlims...)
        if ~(comparison_data isa Nothing) &&
           (Symbol(sn) ∈ propertynames(comparison_data))
            data = getproperty(comparison_data, Symbol(sn))
            data_dates = getproperty(comparison_data, :UTC_datetime)
            spinup_idx = findfirst(spinup_date .<= data_dates)
            CairoMakie.lines!(
                ax,
                data_dates[spinup_idx:end],
                data[spinup_idx:end],
                label = "Data",
                color = "yellow",
            )
        end
        axislegend(ax, position = :lt)
        CairoMakie.save(joinpath(savedir, "$(sn)_$(plot_stem_name).png"), fig)
    end
    return nothing
end

"""
    compare_monthly_fluxes_with_data(
        savedir,
        diagnostics,
        start_date,
        stop_date,
        longlat;
        data_source = "ERA5",
        plot_stem_name = "ERA5",
        spinup_date = start_date,
        short_names = ["lwu", "swu", "shf", "lhf"]
    )

Creates a leaderboard style plot for single column runs where the diagnostics
are stored using the DictWriter, in memory; SHF, LHF, LWu, and SWu are all plotted
against average monthly values from data. This could be extended to pass in other
data sources and compute other flux comparisons.
"""
function compare_monthly_fluxes_with_data(
    savedir,
    diagnostics,
    start_date,
    stop_date,
    longlat;
    data_source = "ERA5",
    plot_stem_name = "ERA5",
    spinup_date = start_date,
    short_names = ["lwu", "swu", "shf", "lhf"],
)
    if data_source == "ERA5"
        data_loader = ERA5DataLoader()
    else
        @error(
            "compare_monthly_fluxes_with_data is currently only configured for ERA5 data."
        )
    end

    source = site_source(diagnostics)
    for sn in short_names
        global_obs_data = get(data_loader, sn)
        ClimaAnalysis.set_reference_date!(global_obs_data, start_date)
        global_obs_data = ClimaAnalysis.window(
            global_obs_data,
            "time",
            left = spinup_date,
            right = stop_date,
        )
        obs_data = ClimaAnalysis.slice(
            global_obs_data,
            lat = longlat[2],
            lon = longlat[1],
        )
        obs_times = obs_data.dims["time"]
        obs_values = obs_data.data
        model_dates, model_output, unit =
            site_timeseries(diagnostics, sn, start_date; source)
        obs_dates = time_to_date.(obs_times, start_date)
        model_spinup_idx = findfirst(spinup_date .<= model_dates)
        obs_spinup_idx = findfirst(spinup_date .<= obs_dates)

        fig = CairoMakie.Figure(size = (800, 400))
        ax = CairoMakie.Axis(
            fig[1, 1],
            xlabel = "Date (UTC)",
            ylabel = "$sn [$(unit)]",
        )
        CairoMakie.lines!(
            ax,
            model_dates[model_spinup_idx:end],
            model_output[model_spinup_idx:end],
            label = "Model",
            color = "blue",
        )
        xlims = extrema(model_dates[model_spinup_idx:end])
        xlims!(ax, xlims...)
        CairoMakie.lines!(
            ax,
            obs_dates[obs_spinup_idx:end],
            obs_values[obs_spinup_idx:end],
            label = "Data",
            color = "yellow",
        )
        axislegend(ax, position = :lt)
        CairoMakie.save(joinpath(savedir, "$(sn)_$(plot_stem_name).png"), fig)
    end
    return nothing
end

"""
    create_partitioning_plots(
        savedir,
        diagnostics,
        start_date,
        stop_date,
        longlat;
        fractions = PARTITION_FRACTIONS,
        era5_to_clima_names = ERA5_PARTITION_TO_CLIMA_NAMES,
    )

Plot the timeseries of each partitioning fraction in `fractions` for a
simulation run at the single site with longitude/latitude pair `longlat`,
saving one png per fraction in `savedir`.

The fractions are the same ones the global partitioning leaderboard shows, so a
site and a global run are described by the same quantities. Every fraction is a
ratio of diagnostics sharing units, so it is formed directly from the recorded
values.

ERA5 is overlaid for the fractions it can supply, loaded under the names
`era5_to_clima_names` maps. Components flagged as prescribed forcing are taken
from the simulation on both sides, since they are the same field there by
construction.

Works with either output writer; see `site_timeseries`.
"""
function create_partitioning_plots(
    savedir,
    diagnostics,
    start_date,
    stop_date,
    longlat;
    fractions = PARTITION_FRACTIONS,
    era5_to_clima_names = ERA5_PARTITION_TO_CLIMA_NAMES,
)
    data_loader = ERA5DataLoader(; era5_to_clima_names)
    available = [d.variable.short_name for d in diagnostics]
    source = site_source(diagnostics)

    for fraction in fractions
        components = union(fraction.numerator, fraction.denominator)
        # A ratio needs every one of its components; there is nothing to plot
        # from a subset of them.
        issubset(components, available) || continue

        # The one place site diagnostics are put on the leaderboard's scale,
        # so nothing downstream has to know whether they still need converting.
        model = Dict(
            c => _site_component(
                site_timeseries(diagnostics, c, start_date; source)...,
                c,
            ) for c in components
        )
        model_dates = model[first(components)].dates
        model_series(names) = sum(model[c].values for c in names)

        fig = CairoMakie.Figure(size = (800, 400))
        ax = CairoMakie.Axis(
            fig[1, 1],
            xlabel = "Date (UTC)",
            ylabel = fraction.ratio_label,
            title = fraction.long_name,
        )
        CairoMakie.lines!(
            ax,
            model_dates,
            model_series(fraction.numerator) ./
            model_series(fraction.denominator),
            label = "Model",
            color = "blue",
        )

        obs_components = setdiff(components, fraction.prescribed)
        # As above: ERA5 supplying part of a ratio is no use, so the overlay is
        # all or nothing.
        if issubset(obs_components, available_vars(data_loader))
            obs = Dict{String, Any}()
            for c in obs_components
                obs[c] = _era5_at_site(data_loader, c, longlat, model_dates)
            end
            # Prescribed components are the same field on both sides, so the
            # simulation supplies them; they are already on the loader's scale.
            for c in fraction.prescribed
                obs[c] = model[c].values
            end
            obs_series(names) = sum(obs[c] for c in names)
            if all(c -> length(obs[c]) == length(model_dates), components)
                CairoMakie.lines!(
                    ax,
                    model_dates,
                    obs_series(fraction.numerator) ./
                    obs_series(fraction.denominator),
                    label = "ERA5",
                    color = "red",
                )
            else
                @info "Skipping the ERA5 overlay for $(fraction.short_name): observations do not cover the simulated months"
            end
        end

        xlims!(ax, extrema(model_dates)...)
        axislegend(ax, position = :lt)
        CairoMakie.save(
            joinpath(savedir, "$(fraction.short_name)_partition.png"),
            fig,
        )
    end
    return nothing
end

"""
    _site_component(dates, values, units, short_name)

Put a site diagnostic on the same footing as the gridded leaderboard: water
fluxes in `mm / day`, the units the ERA5 loader reports them in, and
precipitation as a magnitude. Energy fluxes, and anything already in those
units, pass through untouched.

Returns `(; dates, values, units)` with the units the values ended up in, so
converted output is indistinguishable from output that never needed converting
and a second pass over it does nothing.

This mirrors the `_preprocess_sim_var` methods, which normalize gridded
`OutputVar`s rather than the plain vectors a single site produces.
"""
function _site_component(dates, values, units, short_name)
    # Precipitation is reported downward-positive; the ratios use its magnitude.
    short_name == "precip" && (values = abs.(values))
    units == "kg m^-2 s^-1" &&
        return (; dates, values = values .* 86400.0, units = "mm / day")
    units == "m s^-1" && return (;
        dates,
        values = values .* 1000.0 .* 86400.0,
        units = "mm / day",
    )
    return (; dates, values, units)
end

"""
    _era5_at_site(data_loader, short_name, longlat, dates)

Return the ERA5 timeseries of `short_name` over the months spanned by `dates`,
taken from the grid cell nearest `longlat`: a site sits within one ERA5 cell,
and that cell is what the simulated column is being compared against.
`compare_monthly_fluxes_with_data` reads its observations the same way.
"""
function _era5_at_site(data_loader, short_name, longlat, dates)
    var = get(data_loader, short_name)
    ClimaAnalysis.set_reference_date!(var, first(dates))
    var = ClimaAnalysis.window(
        var,
        "time",
        left = first(dates),
        right = last(dates),
    )
    return ClimaAnalysis.slice(var, lat = longlat[2], lon = longlat[1]).data
end
