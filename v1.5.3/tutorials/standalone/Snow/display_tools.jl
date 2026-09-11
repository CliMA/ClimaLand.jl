using Plots

"""
    cor_map(data, vars)

Visualize the correlation matrix on a DataFrame.

# Arguments
- `data::DataFrame`: the DataFrame of data to investigate.
- `vars::Vector{Symbol}`: the optional list of variables to consider (leaving unspecified visualizes the whole set of variables)
"""
function corr_matrix_plot(
    data::DataFrame;
    vars::Union{Nothing, Vector{Symbol}} = nothing,
)
    df = deepcopy(data)
    if !isnothing(vars)
        df = select(data, vars)
    end
    M = cor(Matrix(df))
    (n, m) = size(M)
    Plots.heatmap(
        M,
        fc = cgrad([:white, :dodgerblue4]),
        xticks = (1:m, names(df)),
        xrot = 90,
        yticks = (1:m, names(df)),
        yflip = true,
    )
    annotate!([
        (j, i, text(round(M[i, j], digits = 3), 8, "Computer Modern", :black))
        for i in 1:n for j in 1:m
    ])
end

"""
    desc(d)

Utility function for checking a DataFrame, equivalent to show(describe(data), allrows = true, allcols = true)

# Arguments
- `d::DataFrame`: the DataFrame of data to describe.
"""
function desc(d::DataFrame)
    show(describe(d), allrows = true, allcols = true)
end

"""
    show_net_weights(params)

Display heatmap of neural network weights by layer.

# Arguments
- `params`: the set of weights to visualize, received by a ModelTools `get_model_ps()` call.
"""
function show_net_weights(params)
    lim = maximum(abs, params)
    display(
        heatmap(
            params,
            aspect_ratio = 1,
            c = cgrad([:red, :white, :blue]),
            clims = (-lim, lim),
        ),
    )
end

"""
    hyperparam_plots(scores, v1, v2)

Display aggregated scores of hyperparameter results, in a heatmap between variables `v1` and `v2`, using the scorign metric/function `sc`.

# Arguments
- `scores::DataFrame`: the set of hyperparameter scores to visualize.
- `v1::Symbol`: The first variable over which to aggregate/visualize scores.
- `v2::Symbol`: The second variable over which to aggregate/visualize scores.
- `sc`: The scoring function, e.g. `mean`, `median`, over which to display scores. Any function taking an AbstractArray and returning a scalar.
"""
function hyperparam_plots(scores::DataFrame, v1::Symbol, v2::Symbol, sc)
    temp = combine(groupby(scores, [v1, v2]), sc => mean, renamecols = false)
    m = unique(temp[!, v1])
    n = unique(temp[!, v2])
    M = zeros((length(m), length(n)))
    for i in 1:length(m)
        for j in 1:length(n)
            M[i, j] =
                temp[(temp[!, v1] .== m[i]) .& (temp[!, v2] .== n[j]), sc][1]
        end
    end
    Plots.heatmap(
        M,
        fc = cgrad([:white, :red]),
        xticks = (1:length(n), unique(temp[!, v2])),
        xlabel = v2,
        xrot = 90,
        yticks = (1:length(m), unique(temp[!, v1])),
        ylabel = v1,
        yflip = true,
    )
    print("Test!")
    annotate!([
        (j, i, text(round(M[i, j], digits = 4), 8, "Computer Modern", :black))
        for i in 1:length(m) for j in 1:length(n)
    ])
end

"""
    siteplot(title, dates, s, l, c; ylab, savename, display_plot)

Display depth timeseries for all series. Utilized for making paper plots. Shades gaps in the data in grey.

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
    #s = [data[!, :z], data[!, :SWE], nn_swes, nn_zs, sn17_zs, sn17_swes, n_zs, sn17o_zs]
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
end

"""
    dplot(title, dates, s, l, c; ylab, yl, savename, display_plot)

Display density timeseries for all series. Utilized for making paper plots. Shades gaps in the data in grey.

# Arguments
- `title::String`: the title to use for the plot.
- `dates::Vector{Union{Date, DateTime}}`: the set of dates from the timeseries generation.
- `s::Vector{Vector{<:Real}}`: the set of series to plot against the dates.
- `l::Vector{Union{Bool, String}}`: the labels to use for each series. a String or "false" is accepted.
- `c::Vector`: the colors to use for each series. Any argument accepted as a color argument for Plots is accepted.
- `ylab::String`: the optional argument to set the y-axis label if depth is not being plotted.
- `yl::Tuple{Real, Real}`: optional argument to set the y-axis scale.
- `savename::String`: the filename used to save the figure if saving is desired. Default is nothing (resulting in no figure saved)
- `display_plot::Bool`: boolean indicating whether to display the plot at the end. Default is `true`.
"""
function dplot(
    title::String,
    dates::Vector{<:Union{Date, DateTime}},
    s::Vector{<:Vector{<:Real}},
    l::Vector{<:Union{Bool, String}} = ["Data", "NN", "SN17", "M", "SN17O"],
    c::Vector = [:black, :dodgerblue, "#33a02c", "#e31a1c", "#6a3d9a"];
    ylab::String = "ρ_snow/ρ_water",
    yl::Tuple{Real, Real} = (0, 1.5),
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
    #s = [tru_d, nn_d, sn17_d, n_d, sn17o_d]
    vspan(skips, color = :grey, alpha = 0.2, label = false)
    for i in 1:length(s)
        vals = Vector{Union{Missing, Real}}(missings(length(ds)))
        temp = allowmissing(deepcopy(s[i]))
        cond = (temp .== 0) .| (isnan.(temp)) .| (isinf.(temp))
        temp[cond] .= missing
        if count(ismissing, temp) == length(temp)
            temp[1] = 0.0
        end
        vals[gd] .= temp
        mc = sum((!).(ismissing.(vals)))
        plabel = (mc > 1) ? l[i] : false
        plot!(ds, vals, color = c[i], label = plabel)
    end
    out = plot!(
        title = title,
        xlabel = "Date",
        ylabel = ylab,
        ylim = yl,
        thickness_scaling = 1.2,
    )
    if display_plot
        display(out)
    end
    if !isnothing(savename)
        savefig(out, savename)
    end
end
