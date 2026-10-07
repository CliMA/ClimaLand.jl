# Compute ClimaLand's FLUXNET2015 RMSE as ILAMB does from the per-site files of
# `fluxnet_ilamb_rmse_sites.jl` and plot it against the ILAMB land-hist models
# of `ilamb_land_hist_site_rmse.csv`. Writes `boxplot_rmse_fluxnet2015.png`,
# `rmse_summary.csv` and `rmse_per_site.csv` to `out/ilamb_rmse/`.
#
# ILAMB's FLUXNET2015 RMSE (`AnalysisMeanStateSites` in ILAMB's `ilamblib.py`)
# is, at each site, the RMSE of the monthly-mean series over the months with
# both model and observations, then the mean over sites. Every model here is
# averaged over the same sites: those where ClimaLand and all land-hist models
# have a value.
#
# Usage:
#   julia --project=.buildkite experiments/integrated/generic_site/fluxnet_ilamb_rmse_plot.jl

import CairoMakie
import DelimitedFiles
import Statistics

const OUTDIR = joinpath(@__DIR__, "out", "ilamb_rmse")
const PANELS = (
    (variable = "hfls", title = "LE"),
    (variable = "hfss", title = "H"),
    (variable = "rsus", title = "SWup"),
    (variable = "rlus", title = "LWup"),
)

"""
    climaland_site_rmse(sites_dir)

Return a `Dict` from `(site, variable)` to ClimaLand's RMSE of the monthly
means at that site, and the number of sites that failed.
"""
function climaland_site_rmse(sites_dir)
    rmse = Dict{Tuple{String, String}, Float64}()
    files = readdir(sites_dir)
    for file in filter(endswith(".csv"), files)
        site = first(splitext(file))
        data, header = DelimitedFiles.readdlm(
            joinpath(sites_dir, file),
            ',';
            header = true,
        )
        isempty(data) && continue
        for variable in unique(data[:, 1])
            rows = (data[:, 1] .== variable) .& isfinite.(Float64.(data[:, 5]))
            any(rows) || continue
            err = Float64.(data[rows, 4]) .- Float64.(data[rows, 5])
            rmse[(site, variable)] = sqrt(Statistics.mean(err .^ 2))
        end
    end
    return rmse, count(endswith(".failed"), files)
end

"""
    land_hist_site_rmse()

Return the model names and a `Dict` from `(site, variable)` to the vector of
per-site RMSEs of the land-hist models (`NaN` where a model has none).
"""
function land_hist_site_rmse()
    data, header = DelimitedFiles.readdlm(
        joinpath(@__DIR__, "ilamb_land_hist_site_rmse.csv"),
        ',',
        String;
        header = true,
    )
    models = String.(vec(header)[3:end])
    values = Dict(
        (data[i, 1], data[i, 2]) =>
            [isempty(v) ? NaN : parse(Float64, v) for v in data[i, 3:end]]
        for i in axes(data, 1)
    )
    return models, values
end

"""
    panel_values(variable, climaland, land_hist)

Return the mean over the common sites of ClimaLand's RMSE and of each
land-hist model's RMSE, and the common sites. Models without any value for
`variable` (JSBACH for rlus) are left out.
"""
function panel_values(variable, climaland, land_hist)
    reported = reduce(
        .|,
        [.!isnan.(v) for ((_, var), v) in land_hist if var == variable],
    )
    sites = sort([
        site for ((site, var), _) in climaland if var == variable &&
            haskey(land_hist, (site, variable)) &&
            all(!isnan, land_hist[(site, variable)][reported])
    ])
    isempty(sites) && return nothing
    others = Statistics.mean(land_hist[(s, variable)][reported] for s in sites)
    clima = Statistics.mean(climaland[(s, variable)] for s in sites)
    return (; others, clima, sites, reported)
end

function bstats(vals)
    q1, med, q3 = Statistics.quantile(vals, (0.25, 0.5, 0.75))
    iqr = q3 - q1
    return (;
        q1,
        med,
        q3,
        lo = max(minimum(vals), q1 - 1.5 * iqr),
        hi = min(maximum(vals), q3 + 1.5 * iqr),
    )
end

# Box of the land-hist models, their jittered values and ClimaLand's value,
# in the style of the global `boxplot_rmse.png` of the leaderboard.
function draw_panel!(fig, col, title, v, y_max, ylabel)
    ax = CairoMakie.Axis(
        fig[1, col],
        title = title,
        titlesize = 18,
        ylabel = ylabel,
        xlabel = "n = $(length(v.sites)) sites\nOthers: $(length(v.others)) models",
        xlabelsize = 12,
        xlabelcolor = :gray,
        xticks = (Float64[], String[]),
        xgridvisible = false,
        limits = (0, 1, 0, y_max),
    )
    CairoMakie.hidexdecorations!(
        ax;
        ticks = true,
        ticklabels = true,
        label = false,
        grid = true,
    )
    cx = 0.40
    box_half = 0.085
    s = bstats(v.others)
    CairoMakie.poly!(
        ax,
        CairoMakie.Rect(cx - box_half, s.q1, 2 * box_half, s.q3 - s.q1);
        color = (:dodgerblue, 0.10),
        strokecolor = :dodgerblue,
        strokewidth = 1.5,
    )
    CairoMakie.lines!(
        ax,
        [cx - box_half, cx + box_half],
        [s.med, s.med];
        color = :darkorange,
        linewidth = 3,
    )
    for (y0, y1) in ((s.lo, s.q1), (s.q3, s.hi))
        CairoMakie.lines!(
            ax,
            [cx, cx],
            [y0, y1];
            color = :dodgerblue,
            linewidth = 1.5,
        )
    end
    for y in (s.lo, s.hi)
        CairoMakie.lines!(
            ax,
            [cx - 0.04, cx + 0.04],
            [y, y];
            color = :dodgerblue,
            linewidth = 1.5,
        )
    end
    for (i, y) in enumerate(v.others)
        jx = cx + ((i * 0.618033) % 1) * 0.16 - 0.08
        CairoMakie.scatter!(
            ax,
            [jx],
            [y];
            color = (:dodgerblue, 0.4),
            strokecolor = :white,
            strokewidth = 0.7,
            markersize = 11,
        )
    end
    CairoMakie.scatter!(
        ax,
        [cx + box_half + 0.18],
        [v.clima];
        color = :firebrick,
        strokecolor = :white,
        strokewidth = 2,
        markersize = 16,
    )
    return ax
end

function plot_rmse(panels, n_failed, outpath)
    y_max = maximum(maximum([v.others; v.clima]) for (_, v) in panels) * 1.18
    fig = CairoMakie.Figure(size = (260 * length(panels) + 140, 600))
    for (col, (p, v)) in enumerate(panels)
        draw_panel!(fig, col, p.title, v, y_max, col == 1 ? "RMSE [W m⁻²]" : "")
    end
    CairoMakie.Label(
        fig[0, :],
        "ClimaLand RMSE — FLUXNET2015 sites (ILAMB definition)";
        fontsize = 18,
        font = :bold,
    )
    legend_elems = [
        CairoMakie.MarkerElement(
            color = (:dodgerblue, 0.4),
            marker = :circle,
            markersize = 11,
        ),
        CairoMakie.LineElement(color = :darkorange, linewidth = 3),
        CairoMakie.PolyElement(
            color = (:dodgerblue, 0.10),
            strokecolor = :dodgerblue,
            strokewidth = 1.5,
        ),
        CairoMakie.MarkerElement(
            color = :firebrick,
            marker = :circle,
            markersize = 14,
        ),
    ]
    CairoMakie.Legend(
        fig[2, :],
        legend_elems,
        [
            "Other models (ILAMB land-hist)",
            "Ensemble median",
            "IQR box",
            "ClimaLand",
        ];
        orientation = :horizontal,
        tellheight = true,
        framevisible = false,
        labelsize = 12,
    )
    note =
        "RMSE as in ILAMB: at each site, RMSE of monthly means against " *
        "FLUXNET2015 (LE_F_MDS, H_F_MDS, SW_OUT, LW_OUT), averaged over the " *
        "sites where every model has a value. ClimaLand is forced by the " *
        "tower meteorology, spun up for one year and scored on the next 12 " *
        "months ($n_failed sites failed or were too short). The other models " *
        "are global runs forced by CRUJRA, GSWP3 or Princeton, sampled at " *
        "the sites and scored over each site's full record.\n" *
        "Disclaimer: For internal use only. A proper Model Intercomparison " *
        "Project (MIP) is beyond scope here, so this comparison is not " *
        "strictly fair: forcings differ across models, tuning targets differ, " *
        "and some variables are prescribed rather than predicted (e.g., LAI " *
        "in ClimaLand). Only a MIP would enable a fair comparison."
    CairoMakie.Label(
        fig[3, :],
        note;
        fontsize = 11,
        color = :gray30,
        word_wrap = true,
        justification = :left,
        halign = :left,
        tellheight = true,
        tellwidth = false,
        padding = (10, 10, 6, 0),
    )
    CairoMakie.save(outpath, fig)
    return outpath
end

function main()
    climaland, n_failed = climaland_site_rmse(joinpath(OUTDIR, "sites"))
    models, land_hist = land_hist_site_rmse()
    panels = []
    open(joinpath(OUTDIR, "rmse_summary.csv"), "w") do io
        println(io, join(["variable", "n_sites", "ClimaLand", models...], ","))
        for p in PANELS
            v = panel_values(p.variable, climaland, land_hist)
            isnothing(v) && continue
            push!(panels, p => v)
            others = fill("", length(models))
            others[v.reported] .= string.(round.(v.others; sigdigits = 4))
            println(
                io,
                join(
                    [
                        p.variable,
                        length(v.sites),
                        round(v.clima; sigdigits = 4),
                        others...,
                    ],
                    ",",
                ),
            )
            @info "$(p.title): ClimaLand $(round(v.clima; sigdigits = 3)), land-hist median $(round(Statistics.median(v.others); sigdigits = 3)) over $(length(v.sites)) sites"
        end
    end
    open(joinpath(OUTDIR, "rmse_per_site.csv"), "w") do io
        println(io, "site,variable,ClimaLand")
        for ((site, variable), r) in sort(collect(climaland))
            println(io, join((site, variable, round(r; sigdigits = 4)), ","))
        end
    end
    isempty(panels) && error("No site has both ClimaLand and land-hist RMSEs")
    plot_rmse(
        panels,
        n_failed,
        joinpath(OUTDIR, "boxplot_rmse_fluxnet2015.png"),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
