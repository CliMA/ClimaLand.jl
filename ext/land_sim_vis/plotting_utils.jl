"""
    make_heatmaps(
        savedir,
        diagdir,
        short_names,
        date;
        plot_name = "figures.pdf",
        levels = nothing,
        plot! = viz.heatmap2D_on_globe!,
        mask = viz.oceanmask(),
        plot_kwargs = Dict(
            :mask => ClimaAnalysis.Utils.kwargs(color = :white),
            :plot => ClimaAnalysis.Utils.kwargs(rasterize = true),
        ),
    )

Generates one .pdf file called `plot_name` in the provided savedir,
this .pdf contains global maps for the provided short_names variables
contained in the provided savedir folder (see ClimaAnalysis documentation),
at the date provided. Variables with multiple levels are plotted at each level
in levels; if levels is nothing, the surface is plotted by default.

The plotting defaults are set for global plots with an ocean mask.
"""
function make_heatmaps(
    savedir,
    diagdir,
    short_names,
    date;
    plot_name = "figures.pdf",
    levels = nothing,
    plot! = viz.heatmap2D_on_globe!,
    mask = viz.oceanmask(),
    plot_kwargs = Dict(
        :mask => ClimaAnalysis.Utils.kwargs(color = :white),
        :plot => ClimaAnalysis.Utils.kwargs(rasterize = true),
    ),
)
    simdir = ClimaAnalysis.SimDir(diagdir)
    mktempdir(savedir) do tmpdir
        for short_name in short_names
            var = get(simdir; short_name)
            avail_dates =
                DateTime(var.attributes["start_date"]) .+
                Second.(ClimaAnalysis.times(var))
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
                            tmpdir,
                            "$(short_name)_$(zval)_$(avail_dates[idx]).pdf",
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
                    joinpath(tmpdir, "$(short_name)_$(avail_dates[idx]).pdf"),
                    fig,
                )
            end
        end
        figures = readdir(tmpdir, join = true)
        pdfunite() do unite
            run(Cmd([unite, figures..., joinpath(savedir, plot_name)]))
        end
    end
    return nothing
end

"""
    check_conservation(savedir, diagdir; plot_name = "conservation_figures.pdf")

Generates one .pdf file called `plot_name`" in the provided 
savedir.

The resulting .pdf contains a time series of the global mean 
(area-weighted) energy and water volume error, in units of 
`J/m^2` and `m`. Only continents are included in the global average. 
"""
function check_conservation(
    savedir,
    diagdir;
    plot_name = "conservation_figures.pdf",
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
    names = ["Global mean energy per area", "Global mean water volume per area"]
    errors = [energy_error, water_volume_error]
    typical_value =
        [@sprintf("%1.2le", mean_energy), @sprintf("%1.2le", mean_water_volume)]
    units = ["J/m²", "m³/m²"]
    mktempdir(savedir) do tmpdir
        for i in 1:2
            fig_cycle = CairoMakie.Figure(size = (600, 400))
            ax = Axis(
                fig_cycle[1, 1],
                xlabel = "Years",
                ylabel = "Global Mean Conservation Error $(units[i])",
                title = "$(names[i]), typical value = $(typical_value[i]) $(units[i])",
            )
            CairoMakie.lines!(ax, times ./ 24 ./ 3600 ./ 365, errors[i])
            CairoMakie.save(joinpath(tmpdir, "$(names[i]).pdf"), fig_cycle)
        end
        figures = readdir(tmpdir, join = true)
        pdfunite() do unite
            run(Cmd([unite, figures..., joinpath(savedir, plot_name)]))
        end
    end
    return nothing
end

"""
    make_ocean_masked_annual_timeseries(
        savedir,
        diagdir,
        short_names;
	plot_name = "annual_timeseries.pdf"
    )

Generates one .pdf file called `plot_name` in the provided savedir,
this .pdf contains the timeseries for the global mean of the provided short_names variables
contained in the provided savedir folder (see ClimaAnalysis documentation). 
"""
function make_ocean_masked_annual_timeseries(
    savedir,
    diagdir,
    short_names;
    plot_name = "annual_timeseries.pdf",
)
    simdir = ClimaAnalysis.SimDir(diagdir)
    mktempdir(savedir) do tmpdir
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
                    round(
                        last(ClimaAnalysis.times(var_sliced) / (365 * 86400)),
                    ), # n years
                )
            ]
            fig_seasonal_cycle = CairoMakie.Figure(size = (600, 400))
            ax = Axis(
                fig_seasonal_cycle[1, 1],
                xlabel = "Month of the year",
                ylabel = "$short_name [$(ClimaAnalysis.units(var))]",
                xticks = (
                    1:1:12,
                    [
                        "Jan",
                        "Feb",
                        "Mar",
                        "Apr",
                        "May",
                        "Jun",
                        "Jul",
                        "Aug",
                        "Sep",
                        "Oct",
                        "Nov",
                        "Dec",
                    ],
                ),
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
                vcat(v, fill(NaN, max_len - length(v))) for
                v in var_global_average
            ] # fill with NaN up to 12, if a vector is shorter
            monthly_averages_var = [
                mean(filter(!isnan, collect([v[i] for v in padded_var])))
                for i in 1:max_len
            ] # compute average of each month, ignoring potential NaNs
            CairoMakie.lines!(
                ax,
                monthly_averages_var,
                color = :black,
                linewidth = 3,
            )
            CairoMakie.xlims!(ax, 1, 12)
            CairoMakie.save(
                joinpath(tmpdir, "$(short_name)_global_monthly.pdf"),
                fig_seasonal_cycle,
            )
        end
        figures = readdir(tmpdir, join = true)
        pdfunite() do unite
            run(Cmd([unite, figures..., joinpath(savedir, plot_name)]))
        end
    end
    return nothing
end
