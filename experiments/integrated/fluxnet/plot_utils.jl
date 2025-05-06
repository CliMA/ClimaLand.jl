"""Plotting utilities for the integrated fluxnet site experiments"""

using Interpolations
using CairoMakie

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
    RMSD = sqrt(mean(abs2, model_hh_avg .- data_hh_avg))
    R² = cor(model_hh_avg, data_hh_avg)^2

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
