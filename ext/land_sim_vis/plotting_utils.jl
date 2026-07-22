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
    N = length(ClimaAnalysis.times(energy_per_area))
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
    for (i, t) in enumerate(times)
        # error = nanmean[(X(t) - X(0) - Expected Change in X)]
        energy_error[i] = ClimaAnalysis.weighted_average_lonlat(
            ClimaAnalysis.apply_oceanmask(
                ClimaAnalysis.slice(energy_per_area, time = t),
            ) - energy_0 - ClimaAnalysis.apply_oceanmask(
                ClimaAnalysis.slice(energy_per_area_change, time = t),
            ),
        ).data[1]

        water_volume_error[i] = ClimaAnalysis.weighted_average_lonlat(
            ClimaAnalysis.apply_oceanmask(
                ClimaAnalysis.slice(water_volume_per_area, time = t),
            ) - water_volume_0 - ClimaAnalysis.apply_oceanmask(
                ClimaAnalysis.slice(water_volume_per_area_change, time = t),
            ),
        ).data[1]

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
            ylabel = "Global Mean Conservation Error [$(units[i])]",
            title = "$(titles[i]), typical value = $(typical_value[i]) $(units[i])",
        )
        CairoMakie.lines!(ax, times ./ 24 ./ 3600 ./ 365, errors[i])
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

        # Can only plot monthly data
        # This assumes EveryCalendarDtSchedule is used with a period of 1 month
        occursin("1 Month", ClimaAnalysis.long_name(var)) || continue

        kwarg_z = ClimaAnalysis.has_altitude(var) ? Dict(:z => 1) : Dict() # if has altitude, take first layer
        var_sliced = ClimaAnalysis.slice(var; kwarg_z...)

        # Get the first and last years of the simulation for windowing
        first_date = first(ClimaAnalysis.dates(var_sliced))
        last_date = last(ClimaAnalysis.dates(var_sliced))
        first_year = year(first_date)
        last_year = year(last_date)

        # var_global_average below is a vector of vector, one for each year of simulation, containing monthly global average of var.
        # i represent a year, from 1 to last year
        # for more details on the ClimaAnalysis functions, see ClimaAnalysis docs.
        var_global_average = [
            ClimaAnalysis.weighted_average_lonlat(
                ClimaAnalysis.apply_oceanmask(
                    ClimaAnalysis.window(
                        var_sliced,
                        "time",
                        left = Dates.DateTime(year, 1),
                        right = Dates.DateTime(year, 12),
                    ),
                ),
            ).data for year in range(first_year, last_year)
        ]

        # Add NaNs for months that are not simulated
        first_year_months =
            ClimaAnalysis.window(
                var_sliced,
                "time",
                left = Dates.DateTime(first_year, 1),
                right = Dates.DateTime(first_year, 12),
            ) |>
            ClimaAnalysis.dates .|>
            Dates.monthname
        indices = indexin(Dates.monthname.(1:12), first_year_months)
        var_global_average[1] =
            [isnothing(i) ? NaN : var_global_average[1][i] for i in indices]

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
    short_names = [d.variable.short_name for d in diagnostics] # short_name_X_average e.g.
    diag_names = [d.output_short_name for d in diagnostics] # short_name_X_average e.g.
    diag_units = [d.variable.units for d in diagnostics]
    for i in 1:length(diag_names)
        dn = diag_names[i]
        unit = diag_units[i]
        sn = short_names[i]
        model_time, model_output = ClimaLand.Diagnostics.diagnostic_as_vectors(
            diagnostics[1].output_writer,
            dn,
        )
        save_Δt = model_time[2] - model_time[1] # in seconds since the start_date. if model_time is an Itime, the epoch should be start_date
        model_dates = time_to_date.(model_time, start_date)
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
    short_names = [d.variable.short_name for d in diagnostics] # short_name
    diag_names = [d.output_short_name for d in diagnostics] # short_name_X_average e.g.
    diag_units = [d.variable.units for d in diagnostics]
    for i in 1:length(diag_names)
        dn = diag_names[i]
        unit = diag_units[i]
        sn = short_names[i]
        model_time, model_output = ClimaLand.Diagnostics.diagnostic_as_vectors(
            diagnostics[1].output_writer,
            dn;
            layer,
        )
        save_Δt = model_time[2] - model_time[1] # in seconds
        model_dates = time_to_date.(model_time, start_date)
        spinup_idx = findfirst(spinup_date .<= model_dates)
        hour_of_day, model_diurnal_cycle = compute_diurnal_cycle(
            model_dates[spinup_idx:end],
            model_output[spinup_idx:end],
        )
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
            hour_of_day, data_diurnal_cycle = compute_diurnal_cycle(
                data_dates[spinup_idx:end],
                data[spinup_idx:end],
            )
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


function compare_monthly_fluxes_with_era5(
    savedir,
    diagnostics,
    start_date,
    stop_date,
    longlat;
    plot_stem_name = "ERA5",
    spinup_date = start_date,
)
    short_names = ["lwu", "swu", "shf", "lhf"] # short_name

    data_loader = ERA5DataLoader()
    for i in 1:length(short_names)
        dn = [
            d.output_short_name for
            d in diagnostics if d.variable.short_name == short_names[i]
        ][1]
        unit = [
            d.variable.units for
            d in diagnostics if d.variable.short_name == short_names[i]
        ][1]
        sn = short_names[i]
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
        model_time, model_output = ClimaLand.Diagnostics.diagnostic_as_vectors(
            diagnostics[1].output_writer,
            dn;
        )
        model_dates = time_to_date.(model_time, start_date)
        model_dates = model_dates .- Month(1)
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


function create_partitioning_plots(savedir, diagnostics)
    diag_short_names = [d.variable.short_name for d in diagnostics]
    # Evaporative fraction
    if "shf" ∈ diag_short_names && "lhf" ∈ diag_short_names
        lhf_dn = [
            d.output_short_name for
            d in diagnostics if d.variable.short_name == "lhf"
        ][1]
        shf_dn = [
            d.output_short_name for
            d in diagnostics if d.variable.short_name == "shf"
        ][1]

        model_time, lhf = ClimaLand.Diagnostics.diagnostic_as_vectors(
            diagnostics[1].output_writer,
            lhf_dn;
        )
        _, shf = ClimaLand.Diagnostics.diagnostic_as_vectors(
            diagnostics[1].output_writer,
            shf_dn;
        )
        start_date = model_time[1].epoch
        model_dates = time_to_date.(model_time, start_date)
        model_dates = model_dates .- Month(1)

        fig = CairoMakie.Figure(size = (800, 400))
        ax = CairoMakie.Axis(
            fig[1, 1],
            xlabel = "Date (UTC)",
            ylabel = "LHF/(SHF+LHF)",
        )
        CairoMakie.lines!(
            ax,
            model_dates,
            lhf ./ (lhf .+ .+shf .+ eps(eltype(shf))),
            label = "Model",
            color = "blue",
        )
        xlims = extrema(model_dates)
        xlims!(ax, xlims...)
        axislegend(ax, position = :lt)
        CairoMakie.save(joinpath(savedir, "evap_fraction.png"), fig)
    end

    # ET partitioning
    if "lhf" ∈ diag_short_names && "clhf" ∈ diag_short_names
        lhf_dn = [
            d.output_short_name for
            d in diagnostics if d.variable.short_name == "lhf"
        ][1]
        clhf_dn = [
            d.output_short_name for
            d in diagnostics if d.variable.short_name == "clhf"
        ][1]

        model_time, lhf = ClimaLand.Diagnostics.diagnostic_as_vectors(
            diagnostics[1].output_writer,
            lhf_dn;
        )
        _, clhf = ClimaLand.Diagnostics.diagnostic_as_vectors(
            diagnostics[1].output_writer,
            clhf_dn;
        )
        start_date = model_time[1].epoch
        model_dates = time_to_date.(model_time, start_date)
        model_dates = model_dates .- Month(1)

        fig = CairoMakie.Figure(size = (800, 400))
        ax = CairoMakie.Axis(fig[1, 1], xlabel = "Date (UTC)", ylabel = "T/ET")
        CairoMakie.lines!(
            ax,
            model_dates,
            clhf ./ (lhf .+ eps(eltype(lhf))),
            label = "Model",
            color = "blue",
        )
        xlims = extrema(model_dates)
        xlims!(ax, xlims...)
        axislegend(ax, position = :lt)
        CairoMakie.save(joinpath(savedir, "et_partition.png"), fig)
    end
    # Runoff partitioning
    if "et" ∈ diag_short_names &&
       "precip" ∈ diag_short_names &&
       "sr" ∈ diag_short_names &&
       "ssr" ∈ diag_short_names
        et_dn = [
            d.output_short_name for
            d in diagnostics if d.variable.short_name == "et"
        ][1]
        sr_dn = [
            d.output_short_name for
            d in diagnostics if d.variable.short_name == "sr"
        ][1]
        ssr_dn = [
            d.output_short_name for
            d in diagnostics if d.variable.short_name == "ssr"
        ][1]
        precip_dn = [
            d.output_short_name for
            d in diagnostics if d.variable.short_name == "precip"
        ][1]
        model_time, et = ClimaLand.Diagnostics.diagnostic_as_vectors(
            diagnostics[1].output_writer,
            et_dn;
        )
        _, precip = ClimaLand.Diagnostics.diagnostic_as_vectors(
            diagnostics[1].output_writer,
            precip_dn;
        )
        _, sr = ClimaLand.Diagnostics.diagnostic_as_vectors(
            diagnostics[1].output_writer,
            sr_dn;
        )
        _, ssr = ClimaLand.Diagnostics.diagnostic_as_vectors(
            diagnostics[1].output_writer,
            ssr_dn;
        )
        start_date = model_time[1].epoch
        model_dates = time_to_date.(model_time, start_date)
        model_dates = model_dates .- Month(1)

        fig = CairoMakie.Figure(size = (800, 400))
        ax = CairoMakie.Axis(fig[1, 1], xlabel = "Date (UTC)", ylabel = "Ratio")
        CairoMakie.lines!(
            ax,
            model_dates,
            et ./ abs.(precip .+ eps(eltype(precip))),
            label = "ET/P",
            color = "blue",
        )
        CairoMakie.lines!(
            ax,
            model_dates,
            (ssr .+ sr) .* 1000 ./ abs.(precip .+ eps(eltype(precip))), # convert to mass flux for runoff
            label = "R/P",
            color = "red",
        )
        xlims = extrema(model_dates)
        xlims!(ax, xlims...)
        axislegend(ax, position = :lt)
        CairoMakie.save(joinpath(savedir, "water_partition.png"), fig)

        fig = CairoMakie.Figure(size = (800, 400))
        ax = CairoMakie.Axis(fig[1, 1], xlabel = "Date (UTC)", ylabel = "Ratio")
        CairoMakie.lines!(
            ax,
            model_dates,
            ssr ./ (ssr .+ sr .+ eps(eltype(sr))),
            label = "SSR/(SSR+SR)",
            color = "blue",
        )

        xlims = extrema(model_dates)
        xlims!(ax, xlims...)
        axislegend(ax, position = :lt)
        CairoMakie.save(joinpath(savedir, "runoff_partition.png"), fig)
    end
    return nothing
end


function create_partitioning_plots(
    savedir,
    diagdir,
    short_names;
    plot! = viz.heatmap2D_on_globe!,
    mask = viz.oceanmask(),
    plot_kwargs = Dict(
        :mask => ClimaAnalysis.Utils.kwargs(color = :white),
        :plot => Dict(
            ClimaAnalysis.Utils.kwargs(rasterize = true)...,
            :colorrange => (0, 1),
        ),
    ),
)
    simdir = ClimaAnalysis.SimDir(diagdir)
    if "shf" ∈ short_names && "lhf" ∈ short_names
        shf_var = get(simdir; short_name = "shf")
        lhf_var = get(simdir; short_name = "lhf")
        # Get the first and last years of the simulation for windowing
        first_date = first(ClimaAnalysis.dates(shf_var))
        last_date = last(ClimaAnalysis.dates(shf_var))
        @assert Second((last_date - first_date)) > Second(86400)
        shf_var = ClimaAnalysis.average_time(
            ClimaAnalysis.window(
                shf_var,
                "time",
                left = last_date - Year(1),
                right = last_date - Month(1),
            ),
        )
        lhf_var = ClimaAnalysis.average_time(
            ClimaAnalysis.window(
                lhf_var,
                "time",
                left = last_date - Year(1),
                right = last_date - Month(1),
            ),
        )
        fig = CairoMakie.Figure(size = (600, 400))
        attributes = Dict(
            "units" => "",
            "short_name" => "ef",
            "long_name" => "Evaporative Fraction",
        )
        ratio = ClimaAnalysis.OutputVar(
            attributes,
            lhf_var.dims,
            lhf_var.dim_attributes,
            lhf_var.data ./ (lhf_var.data .+ shf_var.data),
            lhf_var.dim2index,
            lhf_var.index2dim,
        )
        if mask isa Nothing
            plot!(fig, ratio, more_kwargs = plot_kwargs)
        else
            plot!(fig, ratio, more_kwargs = plot_kwargs, mask = mask)
        end

        CairoMakie.save(joinpath(savedir, "evaporative_fraction.png"), fig)
    end


    if "clhf" ∈ short_names && "lhf" ∈ short_names
        clhf_var = get(simdir; short_name = "clhf")
        lhf_var = get(simdir; short_name = "lhf")
        # Get the first and last years of the simulation for windowing
        first_date = first(ClimaAnalysis.dates(lhf_var))
        last_date = last(ClimaAnalysis.dates(lhf_var))
        @assert Second((last_date - first_date)) > Second(86400)
        clhf_var = ClimaAnalysis.average_time(
            ClimaAnalysis.window(
                clhf_var,
                "time",
                left = last_date - Year(1),
                right = last_date - Month(1),
            ),
        )
        lhf_var = ClimaAnalysis.average_time(
            ClimaAnalysis.window(
                lhf_var,
                "time",
                left = last_date - Year(1),
                right = last_date - Month(1),
            ),
        )
        fig = CairoMakie.Figure(size = (600, 400))
        attributes = Dict(
            "units" => "",
            "short_name" => "tf",
            "long_name" => "Transpiration Fraction",
        )
        ratio = ClimaAnalysis.OutputVar(
            attributes,
            lhf_var.dims,
            lhf_var.dim_attributes,
            clhf_var.data ./ (lhf_var.data .+ 0.01),
            lhf_var.dim2index,
            lhf_var.index2dim,
        )
        if mask isa Nothing
            plot!(fig, ratio, more_kwargs = plot_kwargs)
        else
            plot!(fig, ratio, more_kwargs = plot_kwargs, mask = mask)
        end

        CairoMakie.save(joinpath(savedir, "et_partition.png"), fig)
    end

    if "et" ∈ short_names && "precip" ∈ short_names
        et_var = get(simdir; short_name = "et")
        precip_var = get(simdir; short_name = "precip")
        # Get the first and last years of the simulation for windowing
        first_date = first(ClimaAnalysis.dates(et_var))
        last_date = last(ClimaAnalysis.dates(et_var))
        @assert Second((last_date - first_date)) > Second(86400)
        et_var = ClimaAnalysis.average_time(
            ClimaAnalysis.window(
                et_var,
                "time",
                left = last_date - Year(1),
                right = last_date - Month(1),
            ),
        )
        precip_var = ClimaAnalysis.average_time(
            ClimaAnalysis.window(
                precip_var,
                "time",
                left = last_date - Year(1),
                right = last_date - Month(1),
            ),
        )
        fig = CairoMakie.Figure(size = (600, 400))
        attributes = Dict(
            "units" => "",
            "short_name" => "etp",
            "long_name" => "ET Fraction",
        )
        ratio = ClimaAnalysis.OutputVar(
            attributes,
            et_var.dims,
            et_var.dim_attributes,
            et_var.data ./ abs.(precip_var.data),
            et_var.dim2index,
            et_var.index2dim,
        )
        if mask isa Nothing
            plot!(fig, ratio, more_kwargs = plot_kwargs)
        else
            plot!(fig, ratio, more_kwargs = plot_kwargs, mask = mask)
        end

        CairoMakie.save(joinpath(savedir, "water_partition.png"), fig)
    end

    if "sr" ∈ short_names && "ssr" ∈ short_names
        sr_var = get(simdir; short_name = "sr")
        ssr_var = get(simdir; short_name = "ssr")
        # Get the first and last years of the simulation for windowing
        first_date = first(ClimaAnalysis.dates(sr_var))
        last_date = last(ClimaAnalysis.dates(sr_var))
        @assert Second((last_date - first_date)) > Second(86400)
        sr_var = ClimaAnalysis.average_time(
            ClimaAnalysis.window(
                sr_var,
                "time",
                left = last_date - Year(1),
                right = last_date - Month(1),
            ),
        )
        ssr_var = ClimaAnalysis.average_time(
            ClimaAnalysis.window(
                ssr_var,
                "time",
                left = last_date - Year(1),
                right = last_date - Month(1),
            ),
        )
        fig = CairoMakie.Figure(size = (600, 400))
        attributes = Dict(
            "units" => "",
            "short_name" => "rf",
            "long_name" => "Surface Runoff Fraction",
        )
        ratio = ClimaAnalysis.OutputVar(
            attributes,
            sr_var.dims,
            sr_var.dim_attributes,
            sr_var.data ./ (sr_var.data .+ ssr_var.data),
            sr_var.dim2index,
            sr_var.index2dim,
        )
        if mask isa Nothing
            plot!(fig, ratio, more_kwargs = plot_kwargs)
        else
            plot!(fig, ratio, more_kwargs = plot_kwargs, mask = mask)
        end

        CairoMakie.save(joinpath(savedir, "runoff_partition.png"), fig)
    end

    return nothing
end
