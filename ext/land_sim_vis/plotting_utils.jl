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

function make_diurnal_timeseries(savedir, output_writer, short_names;
                                 plot_name, comparison_data, spinup_dt)
    
    diag_names = collect(keys(output_writer.dict))
    mktempdir(savedir) do tmpdir
        for short_name in short_names
            diag_name = diag_names[findfirst(occursin.(short_name, diag_names))]
            diag_units = ;
            model_time, model_output =  ClimaLand.Diagnostics.diagnostic_as_vectors(output_writer, diag_name; layer = 1) # generalize to 2d
            dt_save = model_time[2] - model_time[1]
            spinup_idx = Int(floor(spinup_dt/dt_save)+1)
            num_days = ;
            if comparison_data isa Nothing
                plot_daily_avg(
                    short_name,
                    model_output[spinup_idx:end],
                    dt_save,
                    num_days,
                    diag_units,
                    savedir,
                    "Model",
                )
            else
                plot_avg_comp(
                    short_name,
                    model_output[spinup_idx:end],
                    dt_save,
                    getproperty(comparison_data, short_name),
                    data_dt,
                    num_days,
                    diag_units
                    savedir,
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

function make_timeseries(savedir, output_writer, short_names;
                                 plot_name, comparison_data, spinup_dt)
    
    diag_names = collect(keys(output_writer.dict))
    mktempdir(savedir) do tmpdir
        for short_name in short_names
            diag_name = diag_names[findfirst(occursin.(short_name, diag_names))]
            diag_units = ;
            model_time, model_output =  ClimaLand.Diagnostics.diagnostic_as_vectors(output_writer, diag_name; layer = 1) # generalize to 2d
            dt_save = model_time[2] - model_time[1]
            spinup_idx = Int(floor(spinup_dt/dt_save)+1)
            num_days = ;
            if comparison_data isa Nothing
                plot!(
                    short_name,
                    model_output[spinup_idx:end],
                    dt_save,
                    num_days,
                    diag_units,
                    savedir,
                    "Model",
                )
            else
                plot!(
                    short_name,
                    model_output[spinup_idx:end],
                    dt_save,
                    getproperty(comparison_data, short_name),
                    data_dt,
                    num_days,
                    diag_units
                    savedir,
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
    

"""Plotting utilities for the integrated fluxnet site experiments"""

using Interpolations
using CairoMakie
using StatsBase

S_PER_DAY = 86400 # Number of seconds in a day

"""This function uses interpolation to convert a time series of data at 
any regular time interval to half hourly data. It takes in the data, the 
vector giving the time stamps of the data in seconds, and the number of days 
in the data series. It returns a vector of the data interpolated to half
hourly intervals."""
function interp_to_hh(data::Vector, timeseries::Vector, num_days::Int64)
    # Rescale the data to a half hourly interval by linear interpolation
    interp = LinearInterpolation(timeseries, data)
    half_hourly = (0.5 * 3600):(0.5 * 3600):(24 * num_days * 3600)
    hh_data = [interp[i] for i in half_hourly]
    return hh_data
end

"""This function takes in a vector of half-hourly data and the number of days 
in the data series and returns the vector giving the diurnal average of 
the data over the time period. (For data which is not half-hourly, use 
interp_hh first.)"""
function hh_to_diurnal_avg(hh_series::Vector, num_days::Int64)
    # Reshape the data into a matrix with each column representing a day
    daily_data = reshape(hh_series, 48, num_days)

    # Average each row over all days to get diurnal average
    hh_avgs = mean(daily_data, dims = 2)
    return hh_avgs
end

"""Given a time series of data, the corresponding time axis, and the number 
of days in the data, linearly interpolates the data to a half hourly time 
scale, and returns a vector giving the diurnal average of the data over the time
period."""
function compute_diurnal_avg(data::Vector, timeseries::Vector, num_days::Int64)
    # Rescale the data to a half hourly interval by linear interpolation
    hh_data = interp_to_hh(data, timeseries, num_days)

    # Get the average of the data over each day
    hh_avgs = hh_to_diurnal_avg(hh_data, num_days)
    return hh_avgs
end

"""This function will be used to plot the average diurnal cycle of a single 
variable. Saves the plot to the directory specified by savedir."""
function plot_daily_avg(
    var_name::String,
    data::Vector,
    data_dt::AbstractFloat,
    num_days::Int64,
    unit::String,
    savedir::String,
    label::String = "data",
)
    # Rescale the data and take the average diurnal cycle
    data_hh_avg =
        compute_diurnal_avg(data, [0:data_dt:(num_days * S_PER_DAY);], num_days)

    # Plot the data diurnal cycle
    fig = CairoMakie.Figure(size = (800, 400))
    ax = CairoMakie.Axis(
        fig[1, 1],
        xlabel = "Hour of day",
        ylabel = "$var_name $(unit)",
        title = "$var_name",
    )
    CairoMakie.lines!(
        ax,
        Array(0.5:0.5:24),
        data_hh_avg[:],
        label = label,
        color = "blue",
    )
    axislegend(ax, position = :lt)
    CairoMakie.save(joinpath(savedir, "$(var_name)_avg.png"), fig)
end

"""This function will be used to plot the comparison of the diurnal average of a 
variable between the model and the data. Saves the plot to the directory
specified by savedir."""
function plot_avg_comp(
    var_name::String,
    model::Vector,
    model_dt::AbstractFloat,
    data::Vector,
    data_dt::AbstractFloat,
    num_days::Int64,
    units::String,
    savedir::String,
)
    # Rescale teh data and take the average diurnal cycle
    model_hh_avg = compute_diurnal_avg(
        model,
        [0:model_dt:(num_days * S_PER_DAY);],
        num_days,
    )
    data_hh_avg =
        compute_diurnal_avg(data, [0:data_dt:(num_days * S_PER_DAY);], num_days)

    # Compute the GOF stats
    RMSD = StatsBase.rmsd(model_hh_avg, data_hh_avg)
    R² = StatsBase.cor(model_hh_avg, data_hh_avg)^2

    # Plot the model and data diurnal cycles
    fig = CairoMakie.Figure(size = (800, 400))
    ax = CairoMakie.Axis(
        fig[1, 1],
        xlabel = "Hour of day",
        ylabel = "$var_name $(units)",
        title = "$var_name: RMSD = $(round(RMSD, digits = 2)), R² = $(round(R²[1][1], digits = 2))",
    )
    CairoMakie.lines!(
        ax,
        Array(0.5:0.5:24),
        model_hh_avg[:],
        label = "Model",
        color = "blue",
    )

    CairoMakie.lines!(
        ax,
        Array(0.5:0.5:24),
        data_hh_avg[:],
        label = "Data",
        color = "yellow",
    )
    axislegend(ax, position = :lt)

    CairoMakie.save(joinpath(savedir, "$(var_name)_avg.png"), fig)
end

"""
This function will be used to plot the comparison of the monthly average of a 
variable between the model and the data. Saves the plot to the directory
specified by savedir.
"""
function plot_monthly_avg_comp(
    var_name::String,
    ref_time,
    model::Vector,
    model_times::Vector,
    data::Vector,
    data_times::Vector,
    units::String,
    savedir::String,
)
    model_avg = compute_monthly_avg(model, model_times)
    data_avg = compute_monthly_avg(data, data_times)

    # Plot the model and data monthly cycles
    fig = CairoMakie.Figure(size = (800, 400))
    ax = CairoMakie.Axis(
        fig[1, 1],
        xlabel = "Month of Year",
        ylabel = "$var_name $(units)",
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
                "Dev",
            ],
        ),
    )

    CairoMakie.lines!(
        ax,
        Array(1:1:12),
        model_avg[:],
        label = "Model",
        color = "blue",
    )

    CairoMakie.lines!(
        ax,
        Array(1:1:12),
        data_avg[:],
        label = "Data",
        color = "yellow",
    )
    axislegend(ax, position = :lt)

    CairoMakie.save(joinpath(savedir, "$(var_name)_monthly_avg.png"), fig)
end

"""
    compute_monthly_avg(data::Vector, times::Vector, ref_date::Dates.DateTime)

Computes the average per month of the `data` measured at `times`, where `times` is in units of seconds past the reference date `ref_date`.
"""
function compute_monthly_avg(
    data::Vector,
    times::Vector,
    ref_date::Dates.DateTime,
)
    months = Dates.month.(ref_date .+ Dates.Second.(times))
    # group by month
    unique_months = unique(months)
    output = zeros(length(unique_months))
    for i in unique_months
        mask = months .== i
        output[i] = mean(data[mask])
    end
    return output, months
end
