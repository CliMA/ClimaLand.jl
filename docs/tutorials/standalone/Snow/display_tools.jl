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
