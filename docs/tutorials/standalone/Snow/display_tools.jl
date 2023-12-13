using Plots

"""
    siteplot(dates, pred, truth, title; LR, savename, display_plot)

Wrapper function for generating pre-formtted plot of predicted vs actual seasonal snow timeseries.

# Arguments
- `dates::Vector{Date}`: the time sequence of the timeseries.
- `pred::Vector{<:Real}`: the predicted values of the timeseries.
- `truth::Vector{<:Real}`: the true values of the timeseries.
- `title`::String`: The title to provide to the plot
- `LR::Vector{<:Real}`: optional specification for inclusion of linear regression predictions as well.
- `savefile::String`: the filename used to save the figure if saving is desired. Default is nothing (resulting in no figure saved)
- `display_plot::Bool`: boolean indicating whether to display the plot at the end. Default is false.
"""
function siteplot(
    dates,
    pred::Vector{<:Real},
    truth::Vector{<:Real},
    title::String = "test";
    LR::Union{Nothing, Vector{<:Real}} = nothing,
    savename = nothing,
    display_plot = false,
)
    plot(
        dates,
        pred,
        label = "Neural Model",
        title = title,
        xlabel = "Date",
        ylabel = "z",
        color = :red,
    )
    if !isnothing(LR)
        plot!(dates, LR, label = "Linear Model", color = :green)
    end
    out = plot!(dates, truth, label = "Raw Data", color = :blue)
    if display_plot
        display(out)
    end
    if !isnothing(savename)
        savefig(out, savename)
    end
end

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
- `params`: the set of weights to visualize.
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
