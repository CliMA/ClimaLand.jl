using Plots, DataFrames, Dates

evaluate(model, input) = model(input)

"""
    custom_loss(x, y, model, n1, n2)

Creates a loss function to be used during training specified by two hyperparameters n1, n2, as outlined in the paper.

# Arguments
- `x`: the input values.
- `y`: output values to compare to.
- `model`: the model over which to evaluate the loss function.
- `n1::Int`: the hyperparameter dictating the scaling of mismatch error.
- `n2::Int`: the hyperparameter dictating the weighting of a given mismatch error by the target magnitude.
"""
function custom_loss(
    x::VecOrMat{FT},
    y::Matrix{FT},
    model,
    n1::Real,
    n2::Real,
)::FT where {FT <: AbstractFloat}
    sum(abs.((model(x) .- y)) .^ n1 .* (1 .+ abs.(y) .^ n2)) / length(y)
end


"""
    make_timeseries(model, timeseries, dt; predictvar, timevar, inputvars, dtype, hole_thresh)

Generate a predicted timeseries given forcing data and the timestep present in that data (holes acceptable).

# Arguments
- `model`: The model used for forecasting (can be any model with a defined `evaluate` call).
- `timeseries::DataFrame`: The input data frame used to generate predictions, including a time variable.
- `dt::Period`: The unit timestep present in the dataframe (i.e. daily dataframe, `dt = Day(1)` or `Second(86400)`).
- `predictvar::Symbol`: The variable to predict from the timeseries. Default is `:z`.
- `timevar::Symbol`: The variable giving the time of each forcing. Default is `:date`.
- `inputvars::Vector{Symbol}`: The variables (in order), to extract from the data to use for predictions.
Default is `[:z, :SWE, :rel_hum_avg, :sol_rad_avg, :wind_speed_avg, :air_temp_avg, :dprecipdt_snow]` like the paper.
- `dtype::Type`: The data type required for input to the model. Default is `Float32`.
- `hole_thresh::Int`: The acceptable number of "holes" in the timeseries for the model to skip over. Default is 30.
"""
function make_timeseries(
    model,
    timeseries::DataFrame,
    dt::Period;
    predictvar::Symbol = :z,
    timevar::Symbol = :date,
    inputvars::Vector{Symbol} = [
        :z,
        :SWE,
        :rel_hum_avg,
        :sol_rad_avg,
        :wind_speed_avg,
        :air_temp_avg,
        :dprecipdt_snow,
    ],
    dtype::Type = Float32,
    hole_thresh::Int = 30,
)
    forcings = Matrix{dtype}(select(timeseries, inputvars))
    pred_idx = findfirst(inputvars .== predictvar)
    check_dates =
        (timeseries[2:end, timevar] - timeseries[1:(end - 1), timevar]) ./ dt
    pred_vals = zeros(nrow(timeseries) - 1)
    pred_series = zeros(nrow(timeseries))
    pred_series[1] = timeseries[1, predictvar]
    countresets = 0
    for j in 2:length(pred_series)
        input = forcings[j - 1, :]
        input[pred_idx] = pred_series[j - 1]
        pred = evaluate(model, input)[1]
        pred_vals[j - 1] = pred
        nperiods = check_dates[j - 1]
        new_val = pred_series[j - 1] + nperiods * Dates.value(Second(dt)) * pred
        pred_series[j] =
            (nperiods <= hole_thresh) ? max(0.0, new_val) :
            timeseries[j, predictvar]  #the "max(0.0, new_val)" is only in the case of traversing holes for a non-negative system, and does not generalize
        if nperiods > hole_thresh
            countresets += 1
        end
    end
    return pred_series, pred_vals, countresets
end


"""
    siteplot(title, dates, s, l, c; ylab, savename, display_plot)

Display depth timeseries for all series. Utilized for making plots in the NeuralDepthModel paper. Shades gaps in the data in grey.

# Arguments
- `title::String`: the title to use for the plot
- `dates::Vector{Union{Date, DateTime}}`: the set of dates from the timeseries generation.
- `s::Vector{Vector{<:Real}}`: the set of series to plot against the dates.
- `l::Vector{Union{Bool, String}}`: the labels to use for each series. a String or "false" is accepted.
- `c::Vector`: the colors to use for each series. Any argument accepted as a color argument for Plots is accepted.
- `ylab::String`: the optional argument to set the y-axis label if depth is not being plotted.
- `savename::String`: the filename used to save the figure if saving is desired. Default is nothing (resulting in no figure saved)
- `display_plot::Bool`: boolean indicating whether to display the plot at the end. Default is `true`.
"""
function siteplot(
    title::String,
    dates::Vector{<:Union{Date, DateTime}},
    s::Vector{<:Vector{<:Union{Missing, Real}}},
    l::Vector = [
        "Data",
        false,
        "NN_SWE",
        "NN_Z",
        "SN17_Z",
        "SN17_SWE",
        "M",
        "SN17O",
    ],
    c::Vector = [
        :black,
        :black,
        :dodgerblue,
        :mediumblue,
        "#33a02c",
        "#b2df8a",
        "#e31a1c",
        "#6a3d9a",
    ];
    ylab = "Snow Depth (m)",
    savename = nothing,
    display_plot::Bool = true,
)
    ds = collect(minimum(dates):Day(1):maximum(dates))
    gd = in.(ds, [dates])
    dts = Dates.value.(dates[2:end] .- dates[1:(end - 1)])
    skips = Vector{Date}()
    for idx in findall(dts .> 1)
        push!(skips, dates[idx])
        push!(skips, dates[idx + 1])
    end
    vspan(skips, color = :grey, alpha = 0.2, label = false)
    for i in 1:length(s)
        vals = Vector{Union{Missing, Real}}(missings(length(ds)))
        vals[gd] .= s[i]
        mc = sum((!).(ismissing.(s[i])))
        plabel = (mc > 1) ? l[i] : false
        plot!(ds, vals, color = c[i], label = plabel)
    end
    out = plot!(
        title = title,
        xlabel = "Date",
        ylabel = ylab,
        thickness_scaling = 1.2,
    )
    if display_plot
        display(out)
    end
    if !isnothing(savename)
        savefig(out, savename)
    end
    return nothing
end
